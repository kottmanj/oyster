import tequila as tq
from contextlib import redirect_stdout
import numpy

class Oyster:

    def __init__(self, geometry, basfilename=None, vqe=None, *args, **kwargs):
        geometry = geometry.strip()
        hydrogens = [a[0].lower()=="h" for a in geometry.split("\n")]
        if not all(hydrogens):
            print("seems like not all atoms are hydrogens!")

        self.geometry = geometry
        if basfilename is None:
            self.basfilename="custombas"
        else:
            self.basfilename=basfilename

        if vqe is None:
            self.vqe = HydrogenChainSPA()
        else:
            self.vqe = vqe

    def __call__(self, basis_set_coefficients, *args, **kwargs):
        n = len(basis_set_coefficients)
        with open(self.basfilename, "w") as f:
            with redirect_stdout(f):
                print('"BASIS "ao basis" SPHERICAL PRINT"')
                print("#BASIS SET: ({n}s) -> [{n}s]".format(n=n))
                for c in basis_set_coefficients:
                    print("H    S")
                    print("{} 1.0".format(c))
                print("END")

        mol = tq.Molecule(geometry=self.geometry, basis_set=self.basfilename)

        if hasattr(self.vqe, "lower"):
            return mol.compute_energy(method=self.vqe.lower())
        else:
            result = self.vqe(mol, *args, **kwargs)

        print(result.energy)
        print(mol.compute_energy(method="fci"))
        return result.energy

class HydrogenChainSPA:

    def __init__(self):
        pass
    def __call__(self,mol,  *args, **kwargs):

        if mol.n_electrons%2 != 0:
            raise Exception("n_electrons={}\nUsed VQE Heuristic is hard coded for even numbers of electrons. Sorry :-(".format(mol.n_electrons))
        # hydrogens
        nh = mol.n_electrons
        # orbitals
        no = mol.n_orbitals
        # orbitals per hydrogen
        oph = no//nh
        # orbitals per pair
        opp = oph*2

        print(no, " ", opp)
        edges = [tuple([i+j for j in range(opp)]) for i in range(0,no,opp)]
        print(edges)

        mol = mol.use_native_orbitals()
        U = mol.make_ansatz(name="HCB-SPA", edges=edges, ladder=False)
        U+= mol.make_ansatz(name="HCB-GD", include_reference=False)

        guess = numpy.eye(no)
        for edge in edges:
            i = edge[0]
            j = edge[opp//2]
            guess[i][j] = 1.0
            guess[j][i] = -1.0

        print(guess)
        result = tq.quantumchemistry.optimize_orbitals(circuit=U, molecule=mol, initial_guess=guess.T, silent=True, use_hcb=True)
        current = result.energy
        for _ in range(10):
            resultx = tq.quantumchemistry.optimize_orbitals(circuit=U, molecule=result.molecule,initial_guess="near_zero", silent=True, use_hcb=True)
            print(resultx.energy)
            if resultx.energy < current:
                result=resultx
            if (resultx.energy - current)<1.e-5:
                break
            current = resultx.energy

        print(result.molecule.integral_manager.orbital_coefficients)
        H = result.molecule.make_hardcore_boson_hamiltonian()
        E = tq.ExpectationValue(H=H, U=U)
        result = tq.minimize(E, silent=True)
        print(result.energy)
        return result

test = Oyster(geometry="H 0.0 0.0 0.0\nH 0.0 0.0 1.0\nH 0.0 0.0 2.0\nH 0.0 0.0 3.0")#, vqe="fci")
test(basis_set_coefficients=[1.0,0.5])



