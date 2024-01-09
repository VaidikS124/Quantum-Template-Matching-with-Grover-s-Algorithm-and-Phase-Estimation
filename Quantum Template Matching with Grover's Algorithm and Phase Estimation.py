%matplotlib inline
# Importing standard Qiskit libraries and configuring account
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, execute, Aer, IBMQ
from qiskit.compiler import transpile, assemble
from qiskit.tools.jupyter import *
from qiskit.visualization import *
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

#need to everything in one circuit for it to plot
matplotlib.rc('text', usetex = True)
matplotlib.rc('font', **{'family' : "sans-serif"})
params= {'text.latex.preamble' : [r'\usepackage{amsmath}']}
figscale = 7
figsize=(figscale*3,figscale*2)
width=0.75
fontsize=34
ticksize=22
color='black'
#data and templates are all in binary form and we compare bit by bit
#define the registers. Completely arbitrary.
#n is the number of qubits(data and templates); nanc is the number of ancilla qubit;
#p is the number of precision; nt is the number of qubits we are not comparing, starting from the 0th
n=6
nanc=1
a=n+nanc
#p=5
nt=1
p=5#int(np.round(np.log2(2*np.pi)+0.5*(n-nt)))
print(p)
#to store the right templates
data = QuantumRegister(n,'data')
#the registers used to check, the last qubit is the ancilla
templates=QuantumRegister(a,'templates')
#the classical register to store the results, of the counting process and the final search process seperately
result_counting=ClassicalRegister(p,'result')
result_answer=ClassicalRegister(n,'answer')
#the counting register
counting=QuantumRegister(p,'counting')
#make the counting register in superposition
search_Circuit=QuantumCircuit(data,templates,counting,result_counting)

#data is 0110
search_Circuit.x(data[1])
search_Circuit.x(data[2])

#search qubits are initialised to 0 while the ancilla is initialised to 1
search_Circuit.x(templates[n])
#Hadmard gate is applied to all qubits
search_Circuit.h(templates[0:a])
#counting qubits are initilised to + as well
search_Circuit.h(counting)
search_Circuit.z(counting[0])
#search_Circuit.x(counting[0])
#search_Circuit.x(counting[1])
search_Circuit.draw(output="mpl")
1126259463#unfortunately qiskit does not have ccccx gate so we have to define
#This is good because we dont have to put in number of qubits by hand
from qiskit.circuit.library.standard_gates import XGate,ZGate
def controlx(n):
    cx_gate = XGate().control(n)
    return cx_gate
def controlz(n):
    cz_gate = ZGate().control(n)
    return cz_gate
#templace matching oracle: one template
#data bit match with template bit
#cx: is a classical AND gate, the result is on the second qubit.  ,templates,data
def template_matching_oracle(qc,ndata):
    '''
    inputregister is the superposition of templates
    data is what we are comparing with
    we are comparing bit by bit, and only when all the comparison is positive we should have 1 in the result register
    qc: the quantum cirsuit this function is performing on
    ndata: number of data
    '''
    qubits=qc.qubits

    ''' #this is the one template case
    for qubit in range(ndata):
        #comparing data with templates qubit by qubit, first is the controll bit
        qc.cx(qubit, qubit+ndata)
        #cx gives a 0 if they match. so x-gates are needed
        qc.x(qubit+ndata)'''
    #this is for two templates
    for qubit in range(1,ndata):
        #comparing data with templates qubit by qubit, first is the controll bit
        qc.cx(qubit, qubit+ndata)
        #cx gives a 0 if they match. so x-gates are needed
        qc.x(qubit+ndata)

    '''
    now is the magic: ancilla qubit which is in the |-> state, with a cccc-x, would if and only if
    adopt a - sign when all bits match
    shame here need to manually change how many control with regard to how many qubits in the template register
    '''
    #qcmin=QuantumCircuit(templates)
    #unfortunately we need to hardcode the control gate in
    '''
    cx_gate=controlx(ndata)
    qc.append(cx_gate,[ndata]+[*range(ndata+1,ndata+ndata+1)])'''
    cx_gate=controlx(ndata-1)
    qc.append(cx_gate,[*range(ndata+1,ndata+ndata+1)])

    #reverse the previous processes to restore the template register
    for qubit in range(1,ndata):
        qc.x(qubit+ndata)
        qc.cx(qubit, qubit+ndata)
    #reverse the ancilla bit
    #qc.h(2*ndata)
    #return (qc)
1126259463#multiple template matching oracle

def template_matching_oracle_multi(qc,ndata,num_match):
    '''
    num_match: number of qubit we are not matching, yielding 2**num_match matching templates.
    qc: the quantum cirsuit this function is performing on
    ndata: number of data
    '''
    #this is the case the from top (qubit0), then that number of qubits are not involved in the match
    for qubit in range(nt,ndata):
        #comparing data with templates qubit by qubit, first is the controll bit
        qc.cx(qubit, qubit+ndata)
        #cx gives a 0 if they match. so x-gates are needed
        qc.x(qubit+ndata)

    '''
    now is the magic: ancilla qubit which is in the |-> state, with a cccc-x, would if and only if
    adopt a - sign when all bits match
    '''
    cx_gate=controlx(ndata-nt)
    #the control-x gate only works on the templates being matched
    qc.append(cx_gate,[*range(ndata+nt,ndata+ndata+1)])

    #reverse the previous processes to restore the template register
    for qubit in range(nt, ndata):
        qc.x(qubit+ndata)
        qc.cx(qubit, qubit+ndata)
#test oracle
oracle_circuit=QuantumCircuit(data,templates)
template_matching_oracle_multi(oracle_circuit,n,nt)
oracle_circuit =transpile(oracle_circuit,basis_gates=['u3','cx','id'],optimization_level=3)
oracle_circuit.draw(output='mpl')

1126259463#diffusion operator
def DiffOpe(qc,ndata):
    #nqubit = len(templates)-1
    #qc=QuantumCircuit(templates)
    #the hadmard
    controlqubits=ndata-1
    cz_gate=controlz(controlqubits)
    for qubit in range(ndata):
        qc.h(ndata+qubit)
        qc.x(ndata+qubit)
    #the phase change, again, we need to hardcode it now
    qc.append(cz_gate,[ndata]+[*range(ndata+1,ndata+ndata)])
    #reverse the process to restore the register
    for qubit in range(ndata):
        qc.x(ndata+qubit)
        qc.h(ndata+qubit)

    #return (qc)

#making Grover's algorithm into a gate
def GroverGate(qc,ndata,num_match):
    #template_matching_oracle(qc,ndata)
    template_matching_oracle_multi(qc,ndata,num_match)
    DiffOpe(qc,ndata)
    return qc

GateCircuit=QuantumCircuit(data,templates)
GateCircuit=GroverGate(GateCircuit,n,nt)
#GateCicuit.draw(output="mpl")
#C_Grover = GateCircuit.to_gate().control()
#print (C_Grover.num_qubits)
G_Grover = GateCircuit.to_gate()
G_Grover.label = "GroverGate"
C_Grover = G_Grover.control()
#cgrit.label = "Grover"
#list(range(1,9))

#the controlled Grover gate part in QPE
iteration=1
#looping through all the counting qubits
for countingqubit in np.arange(n+a,n+a+p):
    for i in np.arange(iteration):
        #this loop is to apply grovers 2^p times
        search_Circuit.append(C_Grover, [countingqubit]+[*range(0,n+a)])
    iteration *=2
#search_Circuit.draw(output="mpl")
#search_Circuit =transpile(search_Circuit,basis_gates=['u3','cx','id'],optimization_level=3)
#search_Circuit.draw(output='mpl')
def qft(n):
    """Creates an n-qubit QFT circuit"""
    circuit = QuantumCircuit(n)
    def swap_registers(circuit, n):
        for qubit in range(n//2):
            circuit.swap(qubit, n-qubit-1)
        return circuit
    def qft_rotations(circuit, n):
        """Performs qft on the first n qubits in circuit (without swaps)"""
        if n == 0:
            return circuit
        n -= 1
        circuit.h(n)
        for qubit in range(n):
            circuit.cp(np.pi/2**(n-qubit), qubit, n)
        qft_rotations(circuit, n)

    qft_rotations(circuit, n)
    swap_registers(circuit, n)
    return circuit
#measuring the counting qubits and output k_t
qft_dagger = qft(p).to_gate().inverse()
qft_dagger.label = "QFT†"
# Do inverse QFT on counting qubits
search_Circuit.append(qft_dagger, list(range(n+a,n+a+p)))
# Measure counting qubits
search_Circuit.measure(list(range(n+a,n+a+p)), result_counting[list(range(p))])
#search_Circuit.draw(output="mpl")
shots=2048
emulator = Aer.get_backend('qasm_simulator')
job = execute(search_Circuit, emulator, shots=shots )
hist = job.result().get_counts()

plot_histogram(hist).savefig('result_n_'+str(n)+'_p_'+str(p)+'_nt_'+str(nt)+'_'+'og'+'.png')

def binary_input(data_ins, n=False):
    #Return the binary representation of the input number as a string
    data_outs=[]
    for data_in in data_ins:
        if n:
            data_outs.append(np.binary_repr(int(data_in), width=n))

        else:
            data_outs.append(np.binary_repr(int(data_in)))
    return np.array(data_outs)

labels=np.arange(2**p)
labels=binary_input(labels,p)
#labels = hist.keys()
probs=[]
print(hist.keys())
for key in labels:
    if np.any(np.array([*hist.keys()])==key):
        probs.append(hist[key]/shots)
    else:
        probs.append(0)
    #print(hist[key])
#print(probs)
print(labels)
fig = plt.figure(figsize=figsize)
plt.bar(labels,probs,width=width,color=color)
plt.xticks(rotation='vertical',fontsize=ticksize)
plt.yticks(rotation='horizontal',fontsize=ticksize)
plt.xlim(-1,len(labels))
plt.xlabel(r'$b$',fontsize=fontsize)
plt.ylabel(r'$p(b)$',fontsize=fontsize)
plt.savefig('result_n_'+str(n)+'_p_'+str(p)+'_nt_'+str(nt)+'_new.png', bbox_inches='tight')
plt.show()

for key in hist.keys():
    if np.any(np.array(hist[key]/shots)>=0.1):
        probs.append(hist[key]/shots)
        print(key,hist[key]/shots)

probs=[]
entry=[]
#print(hist.keys())
for key in hist.keys():
    if np.any(np.array(hist[key]/shots)>=0.001):
        entry.append(key)
        probs.append(hist[key]/shots)
fig = plt.figure(figsize=figsize)
plt.bar(entry,probs,width=width,color=color)
plt.xticks(rotation='vertical')
plt.xlabel(r'\textrm{Template}')
plt.ylabel(r'\textrm{Probability}')
plt.savefig('result_n_'+str(n)+'_p_'+str(p)+'_nt_'+str(nt)+'_clean.png')
plt.show()

measured_str = max(hist, key=hist.get)
measured_int = int(measured_str,2)
print("Register Output = %i" % measured_int)

kt=float(2**(p-2)/measured_int-1/2)
if kt<= 0:
    kt=float(1/(4*(1-measured_int/2**p))-1/2)
print (kt)
#kt1=float(1/(4*(1-measured_int/2**p))-1/2)
#print (kt1)

theta_t=np.pi/(2*(2*kt+1))
num_template=np.sin(theta_t)**2*2**n
print (num_template)
# This is the calculated value of how many times Grover's Gate need to be applied
ktt=float(np.pi/4 *1/np.arcsin(np.sqrt((2**nt)/2**n))-0.5)
#ktt=np.floor(ktt)
print (ktt)
search_Circuit_2=QuantumCircuit(data,templates,result_answer)
search_Circuit_2.x(data[1])
search_Circuit_2.x(data[2])
search_Circuit_2.x(templates[n])
#Hadmard gate is applied to all qubits
search_Circuit_2.h(templates[0:a])
#search_Circuit_2.draw(output="mpl")
#search_circuit=QuantumCircuit(templates,data)
rep=0
if np.round(kt)==0.:
    rep=1
else:
    rep=np.round(kt).astype(int)
for i in range(rep):
    GroverGate(search_Circuit_2,n,nt)

#search_Circuit_2.draw(output='mpl')
#measurement
search_Circuit_2.measure(list(range(n,2*n)), result_answer[list(range(n))])
#experiment with the simulator
shots=2048
emulator = Aer.get_backend('qasm_simulator')
job = execute(search_Circuit_2, emulator, shots=shots )
answer = job.result().get_counts()
#labels = hist.keys()
labels=np.arange(2**n)
labels=binary_input(labels,n)
probs=[]
#print(hist.keys())
for key in labels:
    if np.any(np.array([*answer.keys()])==key):
        probs.append(answer[key]/shots)
    else:
        probs.append(0)
#print(probs)
#print(labels)
fig = plt.figure(figsize=figsize)
plt.bar(labels,probs,width=width,color=color)
plt.xlim(-1,len(labels))
plt.xticks(rotation='vertical',fontsize=ticksize)
plt.yticks(rotation='horizontal',fontsize=ticksize)
plt.xlabel(r'\textrm{Template}',fontsize=fontsize)
plt.ylabel(r'\textrm{Probability}',fontsize=fontsize)
plt.savefig('result_n_'+str(n)+'_p_'+str(p)+'_nt_'+str(nt)+'G'+'_new.png', bbox_inches='tight')
plt.show()

probs=[]
entry=[]
#print(hist.keys())
for key in answer.keys():
    if np.any(np.array(answer[key]/shots)>=0.001):
        entry.append(key)
        probs.append(answer[key]/shots)
        print(answer[key]/shots)
fig = plt.figure(figsize=figsize)
plt.bar(entry,probs,width=width,color=color)
plt.xticks(rotation='vertical')
plt.xlabel(r'\textrm{Template}')
plt.ylabel(r'\textrm{Probability}')
plt.savefig('result_n_'+str(n)+'_p_'+str(p)+'_nt_'+str(nt)+'Gr_clean.png')
plt.show()

plot_histogram(answer).savefig('result_n_'+str(n)+'_p_'+'_nt_'+str(nt)+str(p)+'_'+'og'+'Gr.png')
