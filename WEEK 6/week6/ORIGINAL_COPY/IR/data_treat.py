from IR.IRDipoleApprox import IR

supercell = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

myIR = IR(
    PoscarName="POSCAR",
    BornFileName="BORN",
    ForceConstants=False,
    ForceFileName="FORCE_SETS",
    supercell=supercell,
    nac=False,
    masses=[],
)

myIR.write_spectrum("OscillatorStrengths.txt")
myIR.write_gaussiansmearedspectrum("SmearedOscillatorStrengths.txt", sigma=4)
myIR.plot_spectrum("OscillatorStrengths.pdf")
myIR.plot_gaussiansmearedspectrum("SmearedOscillatorStrengths.pdf", sigma=4)
