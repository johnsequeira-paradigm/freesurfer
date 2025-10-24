# surf2vol - Surface-to-Volume Deformation Using Elastic FEM

## Overview

`surf2vol` is a tool for diffusing surface deformations to volumetric space using elastic finite element methods (FEM). It computes a smooth volumetric transformation field that matches the deformation between two corresponding cortical surfaces, enabling accurate morphing of brain MRI volumes.

**Migrated from PETSc to MFEM (v4.7+)** for improved maintainability and cleaner API.

### Key Features

- **Elastic Finite Element Model**: Uses tetrahedral mesh and linear elasticity
- **Multi-surface Support**: Handle multiple surfaces (white matter, pial, etc.)
- **Topology Preservation**: Automatic detection and correction of mesh inversions
- **Flexible Output**: Supports multiple output formats (MGZ, GCAM, mesh files)
- **Configurable Materials**: Adjust Young's modulus and Poisson ratio
- **Multi-step Refinement**: Progressive mesh refinement for better convergence

## Installation

### Prerequisites

1. **MFEM Library** (v4.7 or later)
   ```bash
   git clone https://github.com/mfem/mfem.git --depth 1 --branch v4.7
   cd mfem
   mkdir build && cd build
   cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/mfem -DMFEM_USE_MPI=NO
   make -j$(nproc)
   make install
   ```

2. **FreeSurfer Libraries**
   - Download from: https://surfer.nmr.mgh.harvard.edu/fswiki/rel7downloads
   - Or build from source with required dependencies

3. **Other Dependencies**
   - CMake 3.10+
   - C++11 compatible compiler
   - TetGen (included in FreeSurfer)
   - BLAS/LAPACK

### Building

```bash
cd fem_elastic
mkdir build && cd build
cmake .. -DMFEM_DIR=$HOME/mfem
make surf2vol
sudo make install  # Optional: install to system
```

## Usage

### Basic Usage

```bash
surf2vol --fixed-mri atlas.mgz \
         --moving-mri subject.mgz \
         --fixed-surf atlas_lh.white \
         --moving-surf subject_lh.white \
         --out morphed.mgz
```

### Multi-Surface Example

```bash
surf2vol --fixed-mri atlas.mgz \
         --moving-mri subject.mgz \
         --fixed-surf lh.white \
         --fixed-surf lh.pial \
         --moving-surf lh.white.moved \
         --moving-surf lh.pial.moved \
         --aparc lh.aparc.annot \
         --aparc lh.aparc.annot \
         --out morphed.mgz \
         --out-mesh transform.tm3d
```

### Advanced Configuration

```bash
surf2vol --fixed-mri atlas.mgz \
         --moving-mri subject.mgz \
         --fixed-surf lh.white \
         --moving-surf lh.white.moved \
         --out output.mgz \
         --fem-steps 5 \
         --elt-vol-min 2 \
         --elt-vol-max 21 \
         --poisson 0.3 \
         --young 10 \
         --compress \
         --cache-transform transform.cache
```

## Command-Line Options

### Required Arguments

| Option | Description |
|--------|-------------|
| `--fixed-mri <file>` | Reference/atlas MRI volume (MGZ format) |
| `--moving-mri <file>` | Subject MRI volume to be morphed (MGZ format) |
| `--fixed-surf <file>` | Reference surface (FreeSurfer format, can specify multiple) |
| `--moving-surf <file>` | Moving surface (FreeSurfer format, can specify multiple) |

**Note:** Number of fixed and moving surfaces must match.

### Optional Arguments

#### Input/Output Options

| Option | Description | Default |
|--------|-------------|---------|
| `-o, --out <file>` | Output morphed volume | `out.mgz` |
| `--out-field <file>` | Output displacement field | `out_field.mgz` |
| `--out-mesh <file>` | Output mesh transformation (.tm3d) | - |
| `--out-surf <root>` | Output surface file root | - |
| `--gcam <file>` | Output GCAM format file | - |
| `--out-affine <file>` | Output affine-only morphed volume | - |
| `--aseg <file>` | ASEG segmentation for atlas | - |
| `--aparc <file>` | APARC annotation (can specify multiple) | - |

#### FEM Parameters

| Option | Description | Default |
|--------|-------------|---------|
| `--elt-vol <value>` | Element volume (sets both min and max) | - |
| `--elt-vol-min <value>` | Minimum element volume | 2.0 |
| `--elt-vol-max <value>` | Maximum element volume | 21.0 |
| `--poisson <value>` | Poisson ratio (must be < 0.5) | 0.3 |
| `--young <value>` | Young's modulus | 10.0 |
| `--fem-steps <n>` | Number of refinement steps | 1 |
| `--fem-end-step <n>` | Final step to compute | -1 |
| `--penalty-weight <value>` | Penalty weight for boundary conditions | 1.0 |

#### Other Options

| Option | Description |
|--------|-------------|
| `-h, --help` | Show help message and exit |
| `--surf-subsample <dist>` | Surface subsampling distance |
| `--cache-transform <file>` | Cache/load linear transform |
| `--compress` | Compress transformation at each step |
| `--topology-old` | Use old topology solver |
| `--use-pial-for-surf` | Use pial surface for mesh construction |
| `--dbg-output <prefix>` | Debug output prefix (saves at each iteration) |

## Algorithm Overview

### Processing Pipeline

1. **Linear Registration**
   - Compute optimal affine transformation using Powell's method
   - Aligns moving surface to fixed surface globally

2. **Mesh Generation**
   - Create tetrahedral mesh using TetGen
   - Progressive refinement from coarse to fine
   - Element size controlled by `--elt-vol-min` and `--elt-vol-max`

3. **Boundary Conditions**
   - Natural BCs: Pin corners of bounding box
   - Multi-point constraints (MFC): Surface correspondence points
   - Weighted penalty method for soft constraints

4. **FEM Solver**
   - Linear elasticity with configurable material properties
   - Conjugate Gradient (CG) solver from MFEM
   - Symmetric system with boundary condition enforcement

5. **Topology Check**
   - Detect inverted tetrahedra
   - Automatic topology correction if needed
   - Choice of old or new topology solver

6. **Volume Morphing**
   - Apply displacement field to moving volume
   - Trilinear interpolation for resampling

### Material Properties

The elastic behavior is controlled by two parameters:

- **Young's Modulus (`--young`)**: Stiffness of the material
  - Higher values → stiffer deformation
  - Lower values → more flexible deformation
  - Typical range: 1.0 - 100.0

- **Poisson Ratio (`--poisson`)**: Incompressibility
  - Must be < 0.5 for stability
  - 0.0 = fully compressible
  - 0.5 = incompressible (use 0.49 max)
  - Typical value: 0.3

### Multi-Step Refinement

Using `--fem-steps N` performs progressive mesh refinement:

```
Step N: Coarse mesh (elt-vol-max)
Step N-1: Medium mesh
...
Step 1: Fine mesh (elt-vol-min)
```

Benefits:
- Better convergence
- Avoids local minima
- Smoother deformation fields

## Input File Formats

### MRI Volumes
- **Format**: MGZ, NIFTI (via FreeSurfer)
- **Requirements**: Same dimensions for fixed and moving
- **Orientation**: Must be in the same space

### Surfaces
- **Format**: FreeSurfer surface format
- **Requirements**:
  - Same number of vertices in corresponding fixed/moving surfaces
  - Surfaces must be in voxel coordinates (automatically converted)
  - Supported types: white matter, pial

### APARC Annotations
- **Format**: FreeSurfer annotation format (.annot)
- **Purpose**: Mask out unreliable regions (e.g., corpus callosum, unknown)
- **Usage**: One annotation per surface (must match surface count)

## Output File Formats

### Morphed Volume (.mgz)
The primary output containing the morphed MRI volume.

### Transformation Mesh (.tm3d)
Custom format storing:
- Tetrahedral mesh structure
- Node displacements
- Can be loaded for future morphing operations

### GCAM Format
GCA morph format compatible with FreeSurfer's mri_warp, etc.

### Displacement Field
Vector field showing displacement at each voxel.

## Examples

### Example 1: Single Surface Morphing

```bash
# Morph subject to atlas using white matter surface
surf2vol --fixed-mri $SUBJECTS_DIR/fsaverage/mri/brain.mgz \
         --moving-mri $SUBJECTS_DIR/subject1/mri/brain.mgz \
         --fixed-surf $SUBJECTS_DIR/fsaverage/surf/lh.white \
         --moving-surf $SUBJECTS_DIR/subject1/surf/lh.white \
         --out subject1_to_atlas.mgz \
         --fem-steps 3
```

### Example 2: Multi-Surface with APARC

```bash
# Use both white and pial surfaces with parcellation
surf2vol --fixed-mri atlas.mgz \
         --moving-mri subject.mgz \
         --fixed-surf lh.white \
         --fixed-surf lh.pial \
         --moving-surf lh.white.moved \
         --moving-surf lh.pial.moved \
         --aparc lh.aparc.annot \
         --aparc lh.aparc.annot \
         --out morphed.mgz \
         --out-mesh mesh.tm3d \
         --fem-steps 5 \
         --compress
```

### Example 3: High-Quality Morphing

```bash
# Fine mesh with multiple refinement steps
surf2vol --fixed-mri atlas.mgz \
         --moving-mri subject.mgz \
         --fixed-surf lh.white \
         --moving-surf lh.white.moved \
         --out high_quality.mgz \
         --elt-vol-min 1.0 \
         --elt-vol-max 25.0 \
         --fem-steps 7 \
         --poisson 0.35 \
         --young 50 \
         --compress \
         --cache-transform transform.cache
```

### Example 4: Debug Mode

```bash
# Save intermediate results at each step
surf2vol --fixed-mri atlas.mgz \
         --moving-mri subject.mgz \
         --fixed-surf lh.white \
         --moving-surf lh.white.moved \
         --out output.mgz \
         --fem-steps 5 \
         --dbg-output debug/morph
# Creates: debug/morph_5.tm3d, debug/morph_4.tm3d, etc.
```

## Performance Considerations

### Memory Usage
- Depends on mesh size (number of tetrahedra)
- Typical: 1-4 GB for standard brain
- Fine meshes (elt-vol-min < 1.0) may require 8+ GB

### Computation Time
- Single step: 5-30 minutes (standard brain)
- Multi-step (5 steps): 30-90 minutes
- Factors affecting time:
  - Number of surface vertices
  - Element volume range
  - Number of refinement steps
  - CPU performance

### Optimization Tips

1. **Start with fewer steps**: Use `--fem-steps 1` for testing
2. **Use coarser meshes**: Increase `--elt-vol-min` for faster computation
3. **Cache transforms**: Use `--cache-transform` to avoid re-computing affine registration
4. **Enable compression**: Use `--compress` to reduce memory for multi-step

## Troubleshooting

### Problem: Topology errors

```
orientation_pbs
found defect for face...
```

**Solution**: The mesh has inverted elements. The tool automatically attempts to fix these. If problems persist:
- Reduce mesh refinement (increase `--elt-vol-min`)
- Use `--topology-old` for alternative solver
- Check surface quality

### Problem: Non-convergence

```
WARNING: MFEM CGSolver did not converge!
```

**Solution**:
- Reduce Young's modulus (`--young 5`)
- Increase element size
- Check boundary conditions (surface correspondence)
- Verify surfaces have same vertex count

### Problem: Out of memory

```
std::bad_alloc
```

**Solution**:
- Increase `--elt-vol-min` (coarser mesh)
- Reduce `--fem-steps`
- Use `--compress` option
- Run on machine with more RAM

### Problem: Missing surfaces

```
Error reading surface...
```

**Solution**:
- Check file paths are correct
- Ensure surfaces are in FreeSurfer format
- Verify read permissions

### Problem: Size mismatch

```
size mismatch for surfaces
```

**Solution**:
- Ensure corresponding fixed/moving surfaces have identical vertex count
- Re-run mris_register or surface generation

## Technical Details

### Finite Element Formulation

Linear elasticity with displacement formulation:
```
∇·σ = 0  (equilibrium)
σ = λtr(ε)I + 2με  (constitutive)
ε = ½(∇u + ∇uᵀ)  (strain)
```

where:
- σ = stress tensor
- ε = strain tensor
- u = displacement field
- λ, μ = Lamé parameters (derived from E and ν)

### Boundary Conditions

1. **Natural (Dirichlet)**: Fixed displacements at boundary nodes
2. **MFC (Multi-Freedom Constraints)**: Penalty method for surface correspondence

### MFEM Solver Configuration

```cpp
mfem::CGSolver cg;
cg.SetRelTol(1.0e-9);      // Relative tolerance
cg.SetMaxIter(10000);       // Maximum iterations
cg.SetPrintLevel(2);        // Verbosity
```

## Comparison: PETSc vs MFEM

| Aspect | PETSc | MFEM |
|--------|-------|------|
| API Style | C-style | Modern C++ |
| Documentation | Extensive | Excellent with examples |
| Learning Curve | Steep | Moderate |
| Community | Large | Growing, very active |
| Code Complexity | High | Low |
| Performance | Excellent | Excellent |

## References

### MFEM
- Website: https://mfem.org/
- Documentation: https://docs.mfem.org/
- GitHub: https://github.com/mfem/mfem
- Examples: https://mfem.org/examples/

### FreeSurfer
- Website: https://surfer.nmr.mgh.harvard.edu/
- Surfaces: https://surfer.nmr.mgh.harvard.edu/fswiki/FsAnatomy
- File Formats: https://surfer.nmr.mgh.harvard.edu/fswiki/FileFormats

### Publications

If you use this tool, please cite:

```
FreeSurfer:
Fischl, B. (2012). FreeSurfer. Neuroimage, 62(2), 774-781.

MFEM:
R. Anderson et al., MFEM: A Modular Finite Element Methods Library,
Computers & Mathematics with Applications, 2021.
```

## Support

- **Issues**: Report bugs on the FreeSurfer GitHub
- **Questions**: FreeSurfer mailing list
- **Documentation**: See `MFEM_MIGRATION.md` for migration details

## License

This tool is part of FreeSurfer and follows the FreeSurfer license.

MFEM is licensed under BSD-3-Clause.

---

**Version**: Post-MFEM Migration
**Last Updated**: 2025
**Maintainer**: FreeSurfer Development Team
