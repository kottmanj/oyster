# use like this
```python
from oyster import Oyster

geometry = "H 0.0 0.0 0.0\nH 0.0 0.0 1.0"
pearl = Oyster(geometry, vqe="fci")
energy = pearl([2.0,1.0,0.5])
```
This evaluates H2 with a basis of 3 Gaussians (exponents 2,1,0.5) per atom.  
If you leave out the `vqe="fci"` part, then an actual VQE is run (not super stable).  

# Dependencies
```bash
pip install tequila-basic
pip install pyscf
pip install qulacs
```
Qulacs is only necessary if you are running a VQE

# Restrictions
currently: Only hydrogen atoms
