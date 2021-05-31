Hands on material for magnetism in DFT using QE


There are 5 directories in which you will find various examples illustrating the course.
You "just" have to execute the various scripts that generate a series of input jobs for QE.
 


Fe_bulk
-------
pw.alat.job : compute total energy and magnetization versus the lattice parameter for Fe.bcc.FM (pw.x)
pw.band.job : compute the spin-polarized band structure of Fe.bcc.FM (pw.x + bands.x)
pw.dos.job  : compute the spin-polarized density of states of Fe.bcc.FM(pw.x + projwfc.x)
pw.fsm.job  : compute fixed spin moment calculation at various lattice parameters of Fe.bcc.FM and Fe.fcc.FM (pw.x)


Fe_surf
--------
pw.dos.job  : compute the spin-polarized density of states for the a 5 layers slab of Febcc(001)(pw.x + projwfc.x)
pw.charge-density.job  : compute charge density (for up and dn spins) of 5 layers slab of Febcc(001)(pw.x + pp.x)

Cr_bulk
--------
pw.alat.job : compute total energy and magnetization versus the lattice parameter for Cr.bcc.AF (pw.x)
pw.fm.job   : compute total energy and magnetization starting of Cr.bcc from a FM configuration (pw.x)
pw.cons.job : constrained (penalized) scf calculation where the magnetization of one of the Cr atom is constrained to its bulk value while the other one is constained to a range of diffferent values
pw.test.job : a test case where the initial magnetization is assymetric +M1 -M2.

Cr_trimer
--------
pw.nc.job  : scf non-collinear calculation of a Cr trimer starting from a frustrated AF collinear magnetization (++-)
pw.nc2.job : scf non-collinear calculation of a Cr trimer starting from an in-plane magnetization convergin towards the ground state.


Fe_wire
------

pw.aniso.alat.job:  calculate the magnetic anisotropy energy Etot(M//)-Etot(Mperp) of Fe wire at two interatomic distances.
pw.aniso.theta.job: calculate the curve Etot(M(theta)) for a Fe wire

