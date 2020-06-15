# MarkovDNA
Markov modeling of DNA replication

This attached function seeks to demonstrate the evolutionary advantage of left-right
symmetry breaking in the kinetic parameters of DNA replica construction
on top of a template strand. The competition between the requirement of
low kinetic barrier for monomer induction and the high kinetic barrier
for monomer retention leads to the symmetry breaking. 

There are two time-scales in the problem: The rate of
H-bonding/dissociation between the free monomers and the template strand,
and the rate of covalent bond formation between monomers on the replica
strand. When the monomer supply is abundant and consequently, the
H-bonding rates are high, there is no evolutionary pressure. But when the
supply is scarce, the H-bonding rates are low and the ability to retain
attached monomers becomes the deciding factor for successful replication.
The DNA replication is modeled as a Markov chain below. '0' represents
the absence of a monomer and '1', its presence. The script models the
growth of a 5-nt long template. Covalent bond formation is not explicitly
included in the chain. It is factored in indirectly as the
retention-advantage factor below, calculated through the probability for
the chain to stay in the '11111' state. 

The DNA replica strand construction is assumed to take place
cooperatively, with the neighboring monomers hydrogen-bonded to the
template influencing the rate of monomer attachment. This rate is not
symmetric w.r.t the left and right neighbors, and can be different. This
difference is encoded in the variable 'a' below, which can be varied from
zero (infinitely high kinetic barrier) to infinity (instantaneous
formation of H-bond). The output produced by the script would show that
it is evolutionarily advantageous for the DNA heteropolymer to have its
left-right symmetry broken, with low barriers to the right (left) for
rapid induction of monomers and high barriers to the left (right) for
retention of the monomers attached to the template. The variable names
are nearly self-explanatory. 


For more information on the model, please refer to "H.  Subramanian, 
R. A. Gatenby, 'Evolutionary advantage of directional symmetry breaking 
in self-replicating polymers', Journal of Theoretical Biology, vol. 446,
pp. 128–136 (2018)".
