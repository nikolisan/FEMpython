# FEM for Python
This project is my implementation of FEM for structures as part of my coursework in MSc Civil and Water Engineering at Cardiff University
This is a project under developement.

## Project Structure
The main package, fempy, is divided into four subpackages:
- **[elements](fempy/elements):** This subpackage contains a set of classes for structural elements, like Triangular2D element, Bar2D element, etc.
- **[mesh](fempy/mesh):** This subpackage provides helper classes for automatic mesh generation.
- **[solvers](fempy/solvers):** This is the main subpackage, which contains all the solver classes for many different applications.
- **[examples](fempy/examples):** As the name states, this package provides some small commented scripts to help the user structure his approach.

## Installation
This project is not yet published to PyPI. To install it you have to build it from source.
#### Build from source
1. Clone the repository and cd inside
```bash
git clone https://github.com/nikolisan/FEMpython.git
cd FEMpython
```
2. Install the library
```bash
python3 -m pip install .
```
3. *Optional.* To uninstall the library
```bash
python3 -m pip uninstall fempy
```

## Usage
After installing the library you can import it to your project
```python
import fempy
# OR
from fempy.elements import TriangularElement2D as Element
from fempy.solvers import PlaneTriangular2D as Solver
```
You can find more usage examples in the [examples](fempy/examples. folder.

<br/>

---

## Authors
* [**Nikolaos Andreakos**](https://github.com/nikolisan) - MSc Civil and Water Engineer
    Find me on [Linked In](https://www.linkedin.com/in/nikolaos-andreakos/)