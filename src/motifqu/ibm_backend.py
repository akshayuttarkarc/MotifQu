"""IBM Quantum hardware backend support.

This module provides functions to run Grover circuits on IBM Quantum
hardware via the IBM Quantum API.

Requirements:
- qiskit-ibm-runtime package
- IBM Quantum API token (get from https://quantum.ibm.com/)
"""

import os
from typing import Optional, Tuple, List
import numpy as np

from qiskit import QuantumCircuit, transpile
from qiskit.primitives import BackendSamplerV2

from motifqu.util import log

# Check for IBM Runtime availability
try:
    from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2
    IBM_RUNTIME_AVAILABLE = True
except ImportError:
    IBM_RUNTIME_AVAILABLE = False


def check_ibm_runtime() -> bool:
    """Check if IBM Quantum Runtime is available."""
    return IBM_RUNTIME_AVAILABLE


def get_ibm_service(token: Optional[str] = None, channel: str = "ibm_quantum_platform") -> "QiskitRuntimeService":
    """Get IBM Quantum Runtime service.
    
    Args:
        token: IBM Quantum API token. If None, uses saved credentials or
               IBMQ_TOKEN environment variable.
        channel: Service channel ('ibm_quantum_platform' or 'ibm_cloud')
    
    Returns:
        QiskitRuntimeService instance
    """
    if not IBM_RUNTIME_AVAILABLE:
        raise ImportError(
            "qiskit-ibm-runtime not installed. Run: pip install qiskit-ibm-runtime"
        )
    
    # Try environment variable if token not provided
    if token is None:
        token = os.environ.get("IBMQ_TOKEN")
    
    if token:
        # Save credentials for future use
        try:
            QiskitRuntimeService.save_account(channel=channel, token=token, overwrite=True)
        except Exception:
            pass  # Already saved or other issue
        return QiskitRuntimeService(channel=channel, token=token)
    else:
        # Try to use saved credentials
        return QiskitRuntimeService(channel=channel)


def list_available_backends(service: "QiskitRuntimeService") -> List[str]:
    """List available IBM Quantum backends.
    
    Args:
        service: QiskitRuntimeService instance
    
    Returns:
        List of backend names
    """
    backends = service.backends()
    return [b.name for b in backends]


def get_best_backend(
    service: "QiskitRuntimeService",
    min_qubits: int = 12,
    simulator_ok: bool = True
) -> str:
    """Get the best available backend for the job.
    
    Args:
        service: QiskitRuntimeService instance
        min_qubits: Minimum qubits required
        simulator_ok: Whether simulators are acceptable
    
    Returns:
        Backend name
    """
    backends = service.backends(
        min_num_qubits=min_qubits,
        simulator=simulator_ok,
        operational=True
    )
    
    if not backends:
        raise RuntimeError(f"No backends available with >= {min_qubits} qubits")
    
    # Sort by queue time (least busy first)
    backends = sorted(backends, key=lambda b: b.status().pending_jobs)
    
    return backends[0].name


def run_on_ibm_hardware(
    qc: QuantumCircuit,
    service: "QiskitRuntimeService",
    backend_name: Optional[str] = None,
    shots: int = 4096,
    optimization_level: int = 3,
    output_dir: Optional[str] = None,
) -> Tuple[np.ndarray, dict]:
    """Run a quantum circuit on IBM Quantum hardware.
    
    Args:
        qc: QuantumCircuit to run (should NOT have save_statevector)
        service: QiskitRuntimeService instance
        backend_name: Specific backend to use (None = auto-select)
        shots: Number of measurement shots
        optimization_level: Transpilation optimization level
        output_dir: Optional directory to save job results
    
    Returns:
        Tuple of (probability distribution as array, raw counts dict)
    """
    import time
    import json
    
    if not IBM_RUNTIME_AVAILABLE:
        raise ImportError("qiskit-ibm-runtime not installed")
    
    n_qubits = qc.num_qubits
    
    # Get backend
    if backend_name is None:
        backend_name = get_best_backend(service, min_qubits=n_qubits)
    
    log(f"Using IBM Quantum backend: {backend_name}")
    backend = service.backend(backend_name)
    
    # Get backend calibration data
    log(f"Fetching calibration data for {backend_name}...")
    try:
        properties = backend.properties()
        if properties:
            log(f"  Calibration date: {properties.last_update_date}")
    except Exception as e:
        log(f"  Could not fetch calibration data: {e}")
    
    # Add measurements to circuit (required for hardware)
    qc_measured = qc.copy()
    qc_measured.measure_all()
    
    # Transpile for the specific backend with optimization
    log(f"Transpiling circuit for {backend_name} (optimization_level={optimization_level})...")
    transpiled = transpile(
        qc_measured,
        backend,
        optimization_level=optimization_level
    )
    
    # Print gate counts after transpilation
    gate_counts = transpiled.count_ops()
    total_gates = sum(gate_counts.values())
    log(f"\n=== Transpilation Results ===")
    log(f"  Circuit depth: {transpiled.depth()}")
    log(f"  Total gates: {total_gates}")
    log(f"  Gate breakdown:")
    for gate, count in sorted(gate_counts.items(), key=lambda x: -x[1]):
        log(f"    {gate}: {count}")
    log(f"  Number of qubits: {transpiled.num_qubits}")
    log(f"=============================\n")
    
    # Run using Sampler primitive
    log(f"Submitting job with {shots} shots...")
    sampler = SamplerV2(backend)
    job = sampler.run([transpiled], shots=shots)
    
    job_id = job.job_id()
    log(f"Job ID: {job_id}")
    log("Waiting for job to complete (checking status every 10 seconds)...\n")
    
    # Poll job status every 10 seconds
    while True:
        status = job.status()
        log(f"  [{time.strftime('%H:%M:%S')}] Job status: {status}")
        
        if status in ["DONE", "COMPLETED"]:
            log(f"\n  Job completed successfully!")
            break
        elif status in ["ERROR", "CANCELLED"]:
            raise RuntimeError(f"Job failed with status: {status}")
        
        time.sleep(10)
    
    result = job.result()
    
    # Extract counts
    pub_result = result[0]
    counts = pub_result.data.meas.get_counts()
    
    # Convert counts to probability array
    N = 2 ** n_qubits
    probs = np.zeros(N)
    
    for bitstring, count in counts.items():
        # Reverse bitstring (qiskit uses little-endian)
        idx = int(bitstring[::-1], 2)
        probs[idx] = count / shots
    
    log(f"\nJob completed. Got {len(counts)} unique outcomes.")
    
    # Save results to output directory if specified
    if output_dir:
        import os
        os.makedirs(output_dir, exist_ok=True)
        
        job_info = {
            "job_id": job_id,
            "backend": backend_name,
            "shots": shots,
            "optimization_level": optimization_level,
            "circuit_depth": transpiled.depth(),
            "total_gates": total_gates,
            "gate_counts": gate_counts,
            "num_qubits": transpiled.num_qubits,
            "counts": counts,
        }
        
        job_file = os.path.join(output_dir, "ibm_job_results.json")
        with open(job_file, "w") as f:
            json.dump(job_info, f, indent=2, default=str)
        log(f"Job results saved to: {job_file}")
    
    return probs, counts


def run_grover_ibm(
    qc: QuantumCircuit,
    token: Optional[str] = None,
    backend_name: Optional[str] = None,
    shots: int = 4096,
    optimization_level: int = 3,
    output_dir: Optional[str] = None,
) -> Tuple[np.ndarray, dict]:
    """Convenience function to run Grover circuit on IBM Quantum.
    
    Args:
        qc: Grover circuit (without measurements or save_statevector)
        token: IBM Quantum API token (optional if saved)
        backend_name: Specific backend (None = auto-select)
        shots: Number of measurement shots
        optimization_level: Transpilation optimization level
        output_dir: Optional directory to save job results
    
    Returns:
        Tuple of (probability array, counts dict)
    """
    service = get_ibm_service(token)
    return run_on_ibm_hardware(
        qc, service, backend_name, shots, optimization_level, output_dir
    )
