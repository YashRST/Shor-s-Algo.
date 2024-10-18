#com\thinking\machines
from qiskit import QuantumCircuit, QuantumRegister
from qiskit_aer import Aer, AerSimulator
from qiskit.circuit.library import QFT, PhaseEstimation
from qiskit.providers import backend
from fractions import Fraction
from math import gcd
from qiskit.transpiler.passes import Depth
from qiskit.quantum_info import Operator
import numpy as np
import pandas as pd
import time
import sys
#from qiskits.connect import IBMQConnector



def initialize_qubits(qc, n_control, n_target):
    """Apply Hadamard gates to the control qubits and initialize the target qubit to |1>."""
    for q in range(n_control):
        qc.h(q)
    qc.x(n_control + n_target - 1)  # Initialize the least significant target qubit to |1>

def controlled_modular_exponentiation(qc, a, N, n_control, n_target):
    """Controlled modular exponentiation."""
    for q in range(n_control):
        exponent = 2 ** q
        for i in range(exponent):
            qc.barrier()  # Barrier for clarity
            qc.cx(q, n_control + (exponent + i) % n_target)  # Control on each target qubit

def apply_iqft(qc, nControl):
    """Apply the inverse QFT to the control qubits."""
    #Apply swaps to reverse the order of qubits
    for qubit in range(nControl // 2):
        qc.swap(qubit, nControl - qubit - 1)
    
    # Apply inverse QFT
    for j in range(nControl):
        qc.h(j)
        for k in range(j + 1, nControl):
            qc.cp(-np.pi / float(2 ** (k - j)), k, j)

def shor(N, a):
    """Main Shor's algorithm function."""
    global numQubit 
    global qcDepth
    n_count = N.bit_length()
    n_control = 2 * n_count  # Number of control qubits
    n_target = n_count

    qc = QuantumCircuit(n_control + n_target, n_control)
    numQubit=n_control+n_target
    initialize_qubits(qc, n_control, n_target)
    controlled_modular_exponentiation(qc, a, N, n_control, n_target)
    apply_iqft(qc, n_control)
    qcDepth=qc.depth()
#    for gate in qc.data:
#        qcDepth = max(qcDepth, gate[1] + 1) 
    qc.measure(range(n_control), range(n_control))
    simulator=AerSimulator()
    job=simulator.run(qc,shots=1024)
    counts=job.result().get_counts(qc)

    # Find the most likely outcome
    max_count = max(counts, key=counts.get)
    measured_phase = int(max_count, 2) / (2 ** n_control)
    r = Fraction(measured_phase).limit_denominator(N).denominator

    # Check if r is a valid period
    if r % 2 == 0 and r != 0:
        factor1 = gcd(a ** (r // 2) - 1, N)
        factor2 = gcd(a ** (r // 2) + 1, N)
        if factor1 != 1 and factor2 != 1:
            return factor1, factor2
    return None


def find_factors(N):
    """Retry with different 'a' values if factors are not found."""
    factors=None
    for a in range(2,N):  # Increased retries for better chances
        if(gcd(a,N)==1): 
            factors = shor(N, a)
        if factors and (factors[0]==1 or factors[1]==1): 
            factors=None
        if factors and factors[0]*factors[1]==N: return factors
    return None

# Example usage
startTime=0
endTime=0
executionTime=0
i=1
print("\tInput\t|"+"\tInput Size\t|"+"\tOutput\t|"+"\tNumber of qubits\t|"+"\tCircuit depth\t|"+"\tExecution time (ms)\t|")
while(i<5):
    primeNumbers=[3,5,7,11,13,17]
    N = primeNumbers[i-1]*primeNumbers[i]  # Try different N values (make sure N is a composite number)
    i+=1
    startTime=time.time_ns()
    factors = find_factors(N)
    endTime=time.time_ns()
    executionTime=(endTime-startTime)/1000000
    
    print("\t "+str(N)+"\t|\t "+str(sys.getsizeof(N))+"\t\t|\t "+str(factors)+"\t|\t\t "+str(numQubit)+"\t\t|\t\t "+str(qcDepth)+"\t|\t "+str(executionTime)+"\t\t|\n")