# SNR-template-attack
We performed a power analysis, we implemented the code for the 3 types of MODELS POWER which are: HD * HW, HW, HD:
HD (Hamming Distance):
There is power consumption when there is a transition in bits i.e. when a change is made from 0 to 1 or vice versa. This expresses the number of FFs that change their value at the next input.
HW (Hamming Weight):
There is power consumption by the number of bits that are 1 regardless of what was before.
HD * HW:
Power consumption is only when there is a bit transition between 0 and 1 i.e. only when the change 0->1 is made.
We only operated on one part of the encryption function, which is the SBOX. The results we get in the CPA section do not tell us how much information from the secret key was leaked from the system, what we get from the CPAs are candidates for the right key and what the probability is the right key among the candidates.
In order to see how much information has been leaked from the system we will use what called MUTUAL INFORMATION so we can know how much information has been leaked in bits.
The advantage of TEMPLATE ATTACK vs. CPA is where TA we need fewer TRACES in order to extract the key right.
