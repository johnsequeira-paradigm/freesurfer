
// STL includes
#include <math.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <iterator>
#include <getopt.h>

// MFEM includes
#include "mfem.hpp"

// FEM includes
#include "solver.h"
#include "fem_3d.h"

// OTHER includes
#include "surf_powell.h" // use Powell to compute the best lin. transform to
// do rigid registration of the surfaces.

// the isfinite defined in utils conflicts with vnl
#undef isfinite

#include "timer.h"
#include "transformUtils.h"
#include "morph.h"

#include "surf_utils.h"
#include "surf_energy.h"

// old topology solver
#include "untangler.h"
// new topology solver
#include "pbCluster_mesh_crop.h"

#include "misc_maths.h"


#define MORPH_WITH_MESH 1
#define DO_MORPH_REGRESSION 0



#include "error.h"
#include "mri.h"
#include "mrinorm.h"
#include "mrisurf.h"
;


const char *Progname;
static char help[] = "Diffuse surface deformation to volumes";

#define USE_SURF_SUBSAMPLING 1


/*--------------------------------------------

Typedefs

------------------------------------------=*/

typedef TSolver<Constructor,3> tSolver;


/*------------------------------------------

CMD-LINE

-------------------------------------------*/

struct IoParams
{
  typedef std::vector<std::string> StringVectorType;

  // fixed volume data
  StringVectorType vstrFixedSurf;
  StringVectorType vstrAparc;
  bool             hasAparc;

  std::string      strFixedMri;
  std::string      strAseg;

  // moving data
  StringVectorType vstrMovingSurf;

  std::string      strMovingMri;


  // output
  std::string      strOutput;
  std::string      strOutputField;
  std::string      strOutputMesh;
  std::string      strOutputSurf; // only root of the name (<index>.surf will be appended)
  std::string      strOutputSurfAffine;
  std::string      strGcam;
  std::string      strOutputAffine;

  // dbg
  std::string      dbgOutput; // will write a morph file at every iteration

  // other options
  //tDblCoords       pixelsPerElt;
  double           eltVolMin, eltVolMax;
  float            poissonRatio;
  float            YoungModulus; 

  bool             compress;

  std::string      strTransform;

  int              iSteps;
  int              iEndStep;

  double           surfSubsample;

  std::string      strDebug;
  bool             bUseOldTopologySolver;
  bool             bUsePialForSurf;

  IoParams();
  int parse(int argc, char* argv[], std::string& errMsg);

  void help_exit();
};




class FunctorTetOrientation
{
public:
  typedef TMesh3d::tElement ElementType;

  FunctorTetOrientation(const ElementType* cpelt) : m_cpelt(cpelt)
  {}

  bool operator()(double dAlpha) const
  {
    return m_cpelt->orientation_test(dAlpha);
  }

  const ElementType* m_cpelt;
};



// Static functions declarations
// !!! initialize PETSC before calling this fct
static
int
process_input_data( const IoParams& params,
                    MRI* &mri_fixed, SurfaceVectorType& vmris_fixed,
                    MRI* &mri_moving, SurfaceVectorType& vmris_moving,
                    MRI* &mri_aseg );


static void
apply_lin_transform(SurfaceVectorType& vmris, float* transform);


static int do_vol_deformation( PointsContainerType& vmris,
                               tSolver& solver,
                               tDblCoords cmin, tDblCoords cmax,
                               int step);

static
float* powell_lin(const SurfaceVectorType& mris_x,
                  const SurfaceVectorType& mris_fx);

static int create_bc_container(PointsContainerType& container,
                               SurfaceVectorType& vmris_fixed,
                               SurfaceVectorType& vmris_moving);

static void compute_fem_error( Transform3SPointer ptransform,
                               const PointsContainerType& container);


// will only use the white surface, if available
//
// the last parameter subsamples the surface - only effective is > 0
static
void
create_surf_container( PointsListType& srcPoints,
                       SurfaceVectorType& vmris,
                       double dist = -1,
                       bool usePialForSurf = false);

// apply transform to container of surface points
//
// This container is used for feeding input to the tetrahedral mesh
//
//
#if MORPH_WITH_MESH
static
void
update_sources(PointsListType& points,
               const CMesh3d* pmesh);
#else
static
void
update_sources(const PointsListType& initPoints,
               PointsListType& srcPoints,
               const Transform3Spointer ptransform);
#endif


unsigned int
update_bc_container(PointsContainerType& container,
                    Transform3SPointer ptransform);


//==========================================================

void check_surface_defects(const SurfaceVectorType& mris_x,
                           const SurfaceVectorType& mris_fx);

void compute_bc_error( const tSolver& solver );

void output_bc_locations( const tSolver& solver,
                          MRI* mri);

//----------------------------------------------------------------
//
// Main function
//
/*

the two datasets are named

ATLAS (mri + wm orig surf)

SRC  (mri + wm orig surf)

The deformation is computed as

ATLAS ==>> SRC

Consequently, the morphing is done in the ATLAS space
(use the direct image to map the deformation of a precise voxel location and then populate the volume)

  Steps in the process:

1. obtain the LS rigid registration between the 2 surfaces. This is currently done using the Powell optimizer.

dir. transform = SRC->ATLAS space
inv. transform = ATLAS->SRC space

2. transport the SRC surface in the atlas space

3. compute the elastic registration from

ATLAS -> SRC

4. Morph the SRC to the ATLAS

This implies

SRC vol(i,j,k) -> dir. transform -> sample using the direct image of the warping

*/

void
dbg_surf2vol(const PointsContainerType& c,
             MRI* mri);

#undef __FUNCT__
#define __FUNCT__ "main"
int
main(int argc,
     char* argv[])
{

  Timer timer;

  // MFEM does not require explicit initialization like PETSc

  // process cmd-line

  IoParams params;

  std::string errMsg;
  params.parse(argc, argv, errMsg);

  if ( !errMsg.empty() )
  {
    std::cerr << " Error parsing the cmd-line\n"
    << errMsg << std::endl;
    exit(1);
  }

  // process the input data

  SurfaceVectorType vmris_fixed;
  SurfaceVectorType vmris_moving;

  MRI*         mri_fixed  = NULL;
  MRI*         mri_aseg   = NULL;
  MRI*         mri_moving    = NULL;

  // if next fct fails -> it exits

  process_input_data( params,
                      mri_fixed, vmris_fixed,
                      mri_moving, vmris_moving,
                      mri_aseg);


  //-------------------------------
  // check 1-2-1 corr. between surf points
  //
  // use resampled surfaces
  if ( vmris_fixed.size() != vmris_moving.size() )
  {
    std::cerr << " Different surface numbers\n";
    exit(1);
  }
  {
    SurfaceVectorType::const_iterator cit_fixed, cit_moving;
    cit_moving = vmris_moving.begin();
    for ( cit_fixed = vmris_fixed.begin();
          cit_fixed != vmris_fixed.end();
          ++cit_fixed, ++cit_moving)
    {
      if ( cit_fixed->mris->nvertices !=
           cit_moving->mris->nvertices )
      {
        std::cerr << " size mismatch for surfaces "
        << cit_fixed - vmris_fixed.begin() << std::endl;
        exit(1);
      }
    } // next cit_fixed, cit_moving
  }


  convert_surf_to_vox( vmris_fixed, mri_fixed);
  convert_surf_to_vox( vmris_moving, mri_moving);

  // init. KSP
  float* transform;
  if ( params.strTransform.empty() )
  {
    std::cout << " applying Powell to get best linear registration\n";
    transform = powell_lin(vmris_moving, vmris_fixed);
  }
  else
  {
    std::cout << " trying to read transform from "
              << params.strTransform <<std::endl;
    transform = read_transform( params.strTransform.c_str() );
    if (!transform)
    {
      std::cout << " reading transform failed -> applying Powell\n";
      transform = powell_lin(vmris_moving, vmris_fixed);
      std::cout << " writing transform to cache file " << params.strTransform
      << std::endl;
      write_transform(transform, params.strTransform.c_str());
    }
  }
  // also compute inverse transform
  float *inv_t;
  inv_transform(transform, &inv_t);
  std::cout << " done with Powell\n";


  // implicitly commit for the working space
  //
  std::cout << " applying linear transform to MOVING surf pos.\n";
  apply_lin_transform(vmris_moving, transform);
  std::cout << "\t DONE\n";

  check_surface_defects(vmris_fixed, vmris_moving);

  // compute the bounding box of the atlas
  //
  // working space = ATLAS
  tDblCoords cmin, cmax;
  tDblCoords TMPcmin, TMPcmax;
  compute_bbox(vmris_fixed, cmin, cmax);
  compute_bbox(vmris_moving, TMPcmin, TMPcmax);
  std::cout << " direct min = " << cmin << " direct max = "
            << cmax << std::endl;
  std::cout << " direct min TMP = " << TMPcmin << " direct max TMP = "
            << TMPcmax << std::endl;
  cmin = cmin - tDblCoords(10.0);
  TMPcmin = TMPcmin - tDblCoords(10.0);

  cmin = max( cmin - tDblCoords(10.0), tDblCoords(0.0) );
  //cmax = cmax + tDblCoords(10.0);
  cmax = tDblCoords(256.0);
  TMPcmax = TMPcmax + tDblCoords(10.0);
  std::cout << " bounding box elts = " << std::endl
  << " min = " << cmin << std::endl
  << " max = " << cmax << std::endl;

  std::cout << " TMP bounding box elts = " << std::endl
  << " TMP min = " << TMPcmin << std::endl
  << " TMP max = " << TMPcmax << std::endl;

  //cmin = min( cmin, TMPcmin );
  //cmax = max( cmax, TMPcmax );

  std::cout << " bounding box elts = " << std::endl
  << " min = " << cmin << std::endl
  << " max = " << cmax << std::endl;

  // due to speed issues, it is much cheaper to
  //     update the BC container at every iteration
  //     however, it is also practical to be able
  //     to compute the error at the very end of the process
  PointsContainerType container, initContainer;
  PointsListType srcPoints, initPoints;

  // encode the deformation field in the atlas surface file
  create_bc_container( initContainer,
                       vmris_fixed,
                       vmris_moving);
  std::copy( initContainer.begin(),
             initContainer.end(),
             std::back_inserter(container) );

  create_surf_container(initPoints,
                        vmris_fixed,
                        params.surfSubsample,
                        params.bUsePialForSurf);
  std::copy( initPoints.begin(),
             initPoints.end(),
             std::back_inserter( srcPoints ) );
  //dbg_surf2vol( initContainer, mri_fixed );
  {
    unsigned int count = 0;
    for ( SurfaceVectorType::iterator it = vmris_moving.begin();
          it != vmris_moving.end(); ++it, ++count )
    {
      MRISsaveVertexPositions( it->mris, TMP_VERTICES );
      convert_vox_to_surf( it->mris, mri_fixed );
      char buf[120];
      sprintf(buf,
              "%s.dbg_surf.%s",
              it->mris->hemisphere?"rh":"lh",
              it->type==white?"white":"pial" );
      std::cout << " writing surface to file "
      << buf << std::endl;
      MRISwrite( it->mris, buf );
      MRISrestoreVertexPositions( it->mris, TMP_VERTICES);
    }
  }

  tDblCoords dbgCoords, dbgImage;
  {
    PointsListType::const_iterator cit = initPoints.begin();
    for ( unsigned int counter = 0;
          counter < 100; ++counter) ++cit;
    dbgCoords = *cit;
    dbgImage = dbgCoords;
  }


  std::cout << " srcPoints size = " << srcPoints.size() << std::endl;

  printf("baseline %g",timer.minutes()); PrintMemUsage(stdout);
  try
  {
    Transform3SPointer ptransform(new gmp::IdentityTransform3d);

    int noSteps = params.iSteps;
    // allow to do some extra steps to finish converging
    for ( int step = noSteps; step > params.iEndStep; --step )
    {
      std::cout << " ======================\n step = " << step  << "\n===============\n";
      printf("step %d %g ",step,timer.minutes());PrintMemUsage(stdout);
      TSolver<Constructor,3> solver;

      // linearly vary the element volume in the given range
      double deltVol = std::max
                       ( params.eltVolMin,
                         params.eltVolMax
                         - (params.eltVolMax - params.eltVolMin) * (noSteps - step) / (noSteps-1)
                       );
      std::cout << "elt_vol= " << deltVol << std::endl;

      DelaunayMesh dmesh(srcPoints, cmin, cmax,deltVol , params.YoungModulus, params.poissonRatio);

      printf("preget %g ",timer.minutes()); PrintMemUsage(stdout);
      CMesh3d* pmesh = dmesh.get();
      printf("postget %g",timer.minutes()); PrintMemUsage(stdout);

      // print basic information about the mesh
      std::cout << " mesh nodes = " << pmesh->get_no_nodes()
      << " mesh elts = " << pmesh->get_no_elts() << std::endl;

      printf("prebuild %g ",timer.minutes()); PrintMemUsage(stdout);
      pmesh->build_index_src();
      printf("postbuild %g ",timer.minutes()); PrintMemUsage(stdout);


      solver.set_mesh(pmesh);

      do_vol_deformation( container, solver, cmin, cmax, std::max(1,step) );

      compute_bc_error(solver);

      // todo - the next timer shows time that could be saved (most of it)
      // namely at the first iteration the time is what it should always be
      Timer t1;

      update_sources( srcPoints, pmesh );
      std::cout << " c0324 time elapsed updating the sources (seconds) = " << t1.seconds() << std::endl;


      if ( !pmesh->orientation_pb(both) )
        std::cout << " NO ORIENTATION PBs FOUND !!!\n";
      else
      {
        std::cout << " orientation_pbs" << std::endl;

        // collect statistics about the orientation pbs
        {
          typedef std::vector<unsigned int> IndexVectorType;
          IndexVectorType vecIdx;
          std::insert_iterator<IndexVectorType> ii(vecIdx,vecIdx.begin());
          pmesh->orientation_pb(dst, ii);

          std::vector<double> vecTime;

          for ( IndexVectorType::const_iterator cit = vecIdx.begin();
                cit != vecIdx.end(); ++cit )
          {
            // retrieve elt
            FunctorTetOrientation fun( pmesh->get_elt( *cit ) );
            // compute max time
            double dtime = FindLeftZeroCrossing(fun, .0,1.);
            vecTime.push_back( dtime );
          }

          //std::ofstream ofs("stats-orientation.txt", std::ios::app );
          //ofs << step << " " << deltVol << " ";
          //std::copy( vecTime.begin(), vecTime.end(),
          //           std::ostream_iterator<double>( ofs, " " ) );
          //ofs << std::endl;
        }

        // solve the orientation pbs
        std::cout << " solving topology problems\n";

        if ( params.bUseOldTopologySolver )
        {
          std::cout << " using old topology solver\n";
          solve_topology_problems( *pmesh, 3 );
        }
        else
        {
          std::cout << " using new topology solver\n";
          TopologySolver( *pmesh, 3 );
        }

        std::cout << "\t DONE\n";
      }

      if ( params.compress )
      {
        std::cout << " compressing morph\n";
        FemTransform3SPointer pTmpTransform(new gmp::FemTransform3d);
        pTmpTransform->set_mesh( std::shared_ptr<TMesh3d>(pmesh) );
        Transform3SPointer pBufTransform(pTmpTransform->convert_to_delta());

        // set a name for the mesh - useful for tracking
        std::ostringstream os;
        os << "compressed_step_" << step;

        update_bc_container(container,
                            pBufTransform);

        std::cout << "dcds mesh direct image = "
        << pmesh->dir_img(dbgImage) << std::endl
        << " dir img for computed transform = "
        << pTmpTransform->img( dbgImage ) << std::endl
        << " dir img for compressed = "
        << pBufTransform->img( dbgImage ) << std::endl;

        pBufTransform->m_strName = os.str();
        pBufTransform->setInitial( ptransform );
        ptransform = pBufTransform;

        dbgImage = ptransform->img(dbgCoords);
        std::cout << "dcds total image for source point = "
        << dbgImage << std::endl;


        ptransform->print();
      }
      else
      {
        FemTransform3SPointer pBufTransform(new gmp::FemTransform3d);

        // set a name for the mesh - useful for tracking
        std::ostringstream os;
        os << "step_" << step;
        pBufTransform->m_strName = os.str();

        pBufTransform->set_mesh( std::shared_ptr<TMesh3d>(pmesh) );

        update_bc_container(container,
                            pBufTransform);

        pBufTransform->setInitial( ptransform );

        ptransform = pBufTransform;
        std::cout << "dcds total image for source point = "
        << ptransform->img(dbgCoords) << std::endl;
      }

      // if dbg option, save the mesh as a volumetric transform
      if ( !params.strDebug.empty() )
      {
        std::shared_ptr<gmp::AffineTransform3d > paffine( new gmp::AffineTransform3d(inv_t) );

        gmp::VolumeMorph morph;
        morph.set_volGeom_fixed( mri_fixed );
        morph.set_volGeom_moving(mri_moving);
        morph.m_transforms.push_back(ptransform);
        morph.m_transforms.push_back(paffine);
        std::ostringstream os;
        os << params.strDebug << "_" << step << ".tm3d";
        std::cout << " saving dbg morph\n";
        morph.save(os.str().c_str());
      }

      compute_fem_error(ptransform,initContainer);

    } // next step
    printf("loopdone ");PrintMemUsage(stdout);

    {
      std::shared_ptr<gmp::AffineTransform3d> paffine(new gmp::AffineTransform3d(inv_t));

      std::cout << " axe\n";
      gmp::VolumeMorph morph;
      morph.m_template = mri_fixed;
      morph.set_volGeom_fixed( mri_fixed );
      morph.set_volGeom_moving( mri_moving );
      morph.m_transforms.push_back(ptransform);
      morph.m_transforms.push_back(paffine);

      MRI* warped = morph.apply_transforms(mri_moving);

      MRIwrite(warped,
               const_cast<char*>( params.strOutput.c_str())  );

      // 02/07/13: orig location has a bug...needed to write it here...
      if( !params.strGcam.empty() ) {
	GCA_MORPH* gcam = morph.exportGcam(mri_moving, FALSE, 0);
	GCAMwrite(gcam, const_cast<char*>( params.strGcam.c_str()));
      }

    }

    std::cout << "baly\n";
    compute_fem_error(ptransform,
                      initContainer);

    if ( !params.strOutputMesh.empty() )
    {
      std::ostringstream os;
      os << params.strOutputMesh << ".tm3d";
      std::cout << " will write volume morph in file "
                << os.str() << std::endl;

      try {
      gmp::VolumeMorph morph;
      morph.set_volGeom_fixed( mri_fixed );
      morph.set_volGeom_moving( mri_moving );
      
      morph.m_transforms.push_back(ptransform);
      std::shared_ptr<gmp::AffineTransform3d>
	paffine( new gmp::AffineTransform3d(inv_t));
      morph.m_transforms.push_back( paffine );	      
      morph.save(os.str().c_str());
      }  
      catch (const char* msg)
	{
	  std::cerr << " exception caught while trying to write transform \n"
		    << msg << std::endl;
	}      
    } 

    // 02/07/13: there is a bug in this; see above for gcam write 
    if (0) {
      if( !params.strGcam.empty() ) {
	// write the m3z version of the morph
	try {
	  std::shared_ptr<gmp::AffineTransform3d> paffine(new gmp::AffineTransform3d(inv_t));
	  gmp::VolumeMorph morph;
	  morph.set_volGeom_fixed( mri_fixed );
	  morph.set_volGeom_moving(mri_moving);
	  morph.m_transforms.push_back(ptransform);
	  morph.m_transforms.push_back(paffine);
	  morph.m_template = mri_fixed;
	  // LZ: test beg
	  //VOL_GEOM vgLike;
	  //initVolGeom(&vgLike);
	  //getVolGeom(morph.m_template, &vgLike);
	  //MRI* mriOut  = morph.apply_transforms(mri_moving, true, &vgLike);
	  //MRI* mriOut  = morph.apply_transforms(mri_moving); 
	  //MRIwrite(mriOut, "tmpout2.mgz");
	  //MRIfree(&mriOut);
	  // LZ: test end
	  GCA_MORPH* gcam = morph.exportGcam(mri_moving, FALSE, 0);
	  GCAMwrite(gcam, const_cast<char*>( params.strGcam.c_str()));
	}
	catch (const char* msg)
	  {
	    std::cerr << " exception caught while trying to write transform \n"
		      << msg << std::endl;
	  }
      }
    }
    // 02/07/13: there is a bug in this; see above for gcam write 
    
    if( !params.strOutputAffine.empty() ) {
      // write the affine morpehd version of the input 
      try {
	std::shared_ptr<gmp::AffineTransform3d> paffine(new gmp::AffineTransform3d(inv_t));
	
	gmp::VolumeMorph morph;
	morph.m_template = mri_fixed;
	morph.set_volGeom_fixed( mri_fixed );
	morph.set_volGeom_moving( mri_moving );
	morph.m_transforms.push_back(paffine);
	
	MRI* warped = morph.apply_transforms(mri_moving);
	
	MRIwrite(warped,
		 const_cast<char*>( params.strOutputAffine.c_str())  );
      }
      catch (const char* msg)
	{
	  std::cerr << " exception caught while trying to write affine morphed input volume \n"
		    << msg << std::endl;
	}
    }
    
  }
  catch (const char* msg)
    {
      std::cout << " caught const char exception \n"
		<< msg << std::endl;
    }
  catch (const std::exception& e)
    {
      std::cout << " Master-net caught exception\n"
		<< e.what() << std::endl;
    }
  
  // MFEM does not require explicit finalization like PETSc

  std::cout << " process performed in " << timer.minutes() << " minutes\n";
  printf("surf2vol done ");PrintMemUsage(stdout);  
  printf("#VMPC# fsurf2vol VmPeak  %d\n",GetVmPeak());

  return 0;
}

///////////////////////////////////////////////////////////////////////////
//
// END MAIN FUNCTION
//
///////////////////////////////////////////////////////////////////////////


//------------------------------------------------------------------
static
void
apply_lin_transform(SurfaceVectorType& vmris, float* transform)
{
  for ( SurfaceVectorType::iterator it = vmris.begin();
        it != vmris.end();
        ++it )
    apply_lin_transform( it->mris, transform);
}

#if MORPH_WITH_MESH

#undef __FUNCT__
#define __FUNCT__ "do_vol_deformation"
static int do_vol_deformation( PointsContainerType& container,
                               tSolver& solver,
                               tDblCoords cmin,
                               tDblCoords cmax,
                               int step)
{
  tDblCoords pt,delta;

  // pin down volume corners to make sure the matrix is positive definite
  delta.set(.0);
  for (int i=0; i<2; ++i)
    for (int j=0; j<2; ++j)
      for (int k=0; k<2; ++k)
      {
        // make sure the boundaries are inside the sampled volume
        pt(0) = i?cmin(0)+.01:cmax(0)-.01;
        pt(1) = j?cmin(1)+.01:cmax(1)-.01;
        pt(2) = k?cmin(2)+.01:cmax(2)-.01;

        solver.add_bc_natural(pt, delta);
      }

  // the container already contains the up-to-date positions of the points
  double dstep = step;
  for ( PointsContainerType::const_iterator pit = container.begin();
        pit != container.end(); ++pit )
  {
    delta = (pit->second - pit->first) / dstep;
    solver.add_bc_mfc(pit->first, delta);
  }
  solver.set_displayLevel(2);
  solver.solve();

  return 0;
}

#else

#undef __FUNCT__
#define __FUNCT__ "do_vol_deformation"
static int do_vol_deformation( PointsContainerType& container,
                               tSolver& solver,
                               tDblCoords cmin,
                               tDblCoords cmax,
                               const Transform3Spointer ptransform,
                               int step)
{
  tDblCoords pt,delta;

  // pin down volume corners to make sure the matrix is positive definite
  delta.set(.0);
  for (int i=0; i<2; ++i)
    for (int j=0; j<2; ++j)
      for (int k=0; k<2; ++k)
      {
        // make sure the boundaries are inside the sampled volume
        pt(0) = i?cmin(0)+.01:cmax(0)-.01;
        pt(1) = j?cmin(1)+.01:cmax(1)-.01;
        pt(2) = k?cmin(2)+.01:cmax(2)-.01;

        solver.add_bc_natural(pt, delta);
      }

  for ( PointsContainerType::const_iterator pit = container.begin();
        pit != container.end();
        ++pit )
  {
    pt = ptransform->img( pit->first );
    delta = (pit->second - pt) / (double)step;
    if ( !pt.isValid() )
    {
      //std::cerr << " invalid pt\n";
      continue;
    }
    solver.add_bc_mfc(pt, delta);
  }
  solver.set_displayLevel(2);
  solver.solve();

  return 0;
}
#endif

static
float* powell_lin(const SurfaceVectorType& vmris_x,
                  const SurfaceVectorType& vmris_fx)
{
  PointsContainerType container;

  VERTEX* pvtx_x = NULL;
  VERTEX* pvtx_fx = NULL;
  unsigned int nvertices;
  tDblCoords pt_x, pt_fx;

  int linInc = 1; // Default value, no longer using PETSc command-line options

  SurfaceVectorType::const_iterator cit_x, cit_fx;
  cit_fx = vmris_fx.begin();

  for ( cit_x = vmris_x.begin();
        cit_x != vmris_x.end();
        ++cit_x, ++cit_fx )
  {
    MRI_SURFACE* mris_x = cit_x->mris;
    MRI_SURFACE* mris_fx = cit_fx->mris;

    pvtx_x = &( mris_x->vertices[0] );
    pvtx_fx = &( mris_fx->vertices[0] );
    nvertices = (unsigned int)mris_x->nvertices;

    for (unsigned int ui=0;
         ui < nvertices;
         ui+=linInc, pvtx_x+=linInc, pvtx_fx+=linInc )
    {
      // don't use medial or unknown regions in linear fit
      if ( pvtx_x->ripflag ) continue;

      pt_x(0) = pvtx_x->x;
      pt_x(1) = pvtx_x->y;
      pt_x(2) = pvtx_x->z;

      pt_fx(0) = pvtx_fx->x;
      pt_fx(1) = pvtx_fx->y;
      pt_fx(2) = pvtx_fx->z;

      container.push_back( std::make_pair(pt_x, pt_fx));
    } // next ui, pvtx_x, pvtx_fx
  } // next cit_x, cit_fx

  // setup initial transform as identity
  float *transform = new float[12]; // line storage
  std::fill( transform, transform + 12, 0.0f);
  for ( int i=0; i<3; ++i) transform[ 4*i ] = 1.0f;

  // call routine
  std::cout << " energy before calling Powell = "
  << energy( container, transform) << std::endl;
  powell_minimize( container,
                   transform,
                   energy);
  std::cout << " energy after Powell = " << energy( container, transform ) << std::endl;

  return transform;
}

static int create_bc_container(PointsContainerType& container,
                               SurfaceVectorType& vmris_fixed,
                               SurfaceVectorType& vmris_moving)
{
  double dirty = 1.0; // Default value, no longer using PETSc command-line options
  std::cout << " DIRTY value = " << dirty << std::endl;

  tDblCoords pt, img;
  container.clear();

  SurfaceVectorType::const_iterator cit_fixed, cit_moving;
  cit_moving = vmris_moving.begin();
  for ( cit_fixed = vmris_fixed.begin();
        cit_fixed != vmris_fixed.end();
        ++cit_fixed, ++cit_moving)
  {
    MRI_SURFACE* mris_fixed = cit_fixed->mris;
    MRI_SURFACE* mris_moving = cit_moving->mris;
#if 0
    int markedVertices = MRISsubsampleDist(mris_fixed, 2.0f);
#endif
    VERTEX* pvtx_fixed = &( mris_fixed->vertices[0] );
    unsigned int nvertices = (unsigned int)mris_fixed->nvertices;
    VERTEX* pvtx_moving = &( mris_moving->vertices[0]);

    unsigned int ripCounter = 0;
    for ( unsigned int ui = 0;
          ui < nvertices;
          ++ui, ++pvtx_fixed, ++pvtx_moving)
    {
      if ( pvtx_fixed->ripflag )
      {
        ++ripCounter;
        continue;
      }
      // don't subsample the BC container
#if 0
      if ( !pvtx_fixed->val ) continue;
#endif
      pt(0) = pvtx_fixed->x;
      pt(1) = pvtx_fixed->y;
      pt(2) = pvtx_fixed->z;

      img(0) = pvtx_moving->x;
      img(1) = pvtx_moving->y;
      img(2) = pvtx_moving->z;

      container.push_back( std::make_pair(pt, img) );
    } // next ui, pvtx_fixed, pvtx_moving
    std::cout << " create BC container - ripCounter = " << ripCounter << std::endl;
  } // next cit_fixed, cit_moving

  return 0;
}

static
void
create_surf_container( PointsListType& srcPoints,
                       SurfaceVectorType& vmris,
                       double subsampleDist,
                       bool usePialForSurf)
{
  tDblCoords pt;

  std::cout << " creating container for the mesher\n";

  for ( SurfaceVectorType::const_iterator cit = vmris.begin();
        cit != vmris.end(); ++cit )
  {
    if ( cit->type == pial || !usePialForSurf )
      continue;

    MRI_SURFACE* mris = cit->mris;
    // subsample the surface

    if ( subsampleDist > 0 )
    {
      int markedVertices = MRISsubsampleDist(mris, 3.5f);
      std::cout << " create_surf_container - init size = " << mris->nvertices
      << " subsampled = " << markedVertices << std::endl;
    }
    VERTEX* pvtx = &( mris->vertices[0] );
    const unsigned int nvertices = mris->nvertices;
    for (unsigned int ui=0;
         ui < nvertices; ++ui, ++pvtx)
    {
      if ( pvtx->ripflag ) continue; // uses the aparc
      if ( subsampleDist>0 && !pvtx->val ) continue;

      pt(0) = pvtx->x;
      pt(1) = pvtx->y;
      pt(2) = pvtx->z;

      srcPoints.push_back( pt );
    }
  }
}

#if MORPH_WITH_MESH
struct PointMorphWithMesh
{
  PointMorphWithMesh(const CMesh3d* pm=NULL) : pmesh(pm)
  {}
  const CMesh3d* pmesh;
  tDblCoords operator()(tDblCoords pt) const
  {
    return pmesh->dir_img(pt);
  }
};
struct InvalidPointFilter
{
  bool operator()(tDblCoords pt) const
  {
    return !pt.isValid();
  }
};

static
void
update_sources(PointsListType& points,
               const CMesh3d* pmesh)
{
  std::transform( points.begin(), points.end(),
                  points.begin(), PointMorphWithMesh(pmesh) );
  std::remove_if( points.begin(), points.end(), InvalidPointFilter() );

}
#else
static
void
update_sources(const PointsListType& initPoints,
               PointsListType& srcPoints,
               const Transform3Spointer ptransform)
{
  tDblCoords pt;

  srcPoints.clear();

  unsigned int count = 0;
  for ( PointsListType::const_iterator pit = initPoints.begin();
        pit != initPoints.end(); ++pit )
  {
    pt = ptransform->img( *pit );
    if ( !pt.isValid() ) ++count;
    srcPoints.push_back(pt);
  } // next pit

  std::cout << " srcPoints size = " << srcPoints.size() << std::endl;

  if ( count )
    std::cerr << " errors while updating sources -> count = "
    << count << std::endl;
}
#endif

#ifdef USE_COMPUTE_MORPH
#if 1
static MRI* compute_morph( Transform3SPointer ptransform,
                           MATRIX* mat_lt,
                           MRI* mri_x,
                           MRI* mri_atlas,
                           MRI* mri_src,
                           tDblCoords cmin,
                           tDblCoords cmax
                         )
{
  MRI* warped = NULL;
  // recover pointer for the transform
  float* plin = new float[12];
  for ( int i=0; i<4; ++i)
    for ( int j=0; j<3; ++j)
      plin[i*3+j] = mat_lt->rptr[j+1][i+1];

  try
  {
    std::shared_ptr<gmp::Transform<3> > bpTransform(ptransform);
    std::shared_ptr<gmp::AffineTransform3d>
    paffine( new gmp::AffineTransform3d(plin));

    gmp::VolumeMorph volMorph;
    volMorph.m_template = mri_atlas;

    volMorph.m_transforms.push_back( bpTransform );
    volMorph.m_transforms.push_back( paffine );

    warped = volMorph.apply_transforms(mri_src);

    std::cout << " after transform\n";
  }
  catch (const gmpErr& excp)
  {
    std::cerr << "Exception caught -> "
    << excp.what() << std::endl;
    exit(1);
  }


  delete[] plin;

  return warped;

}
#else

static MRI* compute_morph( const CMesh3d& mesh,
                           MATRIX* mat_lt,
                           MRI* mri_x,
                           MRI* mri_atlas,
                           MRI* mri_src,
                           tDblCoords cmin,
                           tDblCoords cmax)
{

  MRI* warped = NULL;
  std::cout << " 1\n";
  // recover pointer for the transform
  float* plin = new float[12];
  for ( int i=0; i<4; ++i)
    for ( int j=0; j<3; ++j)
      plin[i*3+j] = mat_lt->rptr[j+1][i+1];

  std::cout << " 2\n";

  try
  {
    gmp::AffineTransform3d affine(plin);
    gmp::FemTransform3d femTransform(&mesh);
    femTransform.set_bbox( cmin, cmax);

    std::cout << " 3\n";

    gmp::VolumeMorph volMorph;
    volMorph.m_template = mri_atlas;
    volMorph.m_input = mri_src;

    std::cout << " got here\n";
    volMorph.m_transforms.push_back( &femTransform );

    warped = volMorph.apply_transforms();
    std::cout << " after transform\n";
  }
  catch (const gmpErr& excp)
  {
    std::cerr << "Exception caught -> "
    << excp.what() << std::endl;
    exit(1);
  }


  delete[] plin;

  return warped;

}

#endif
#endif //ifdef USE_COMPUTE_MORPH


/*-----------------------------------

CMD-LINE

-----------------------------------*/


IoParams::IoParams()
    : vstrFixedSurf(),
    vstrAparc(),
    hasAparc(false),
    strFixedMri(),
    strAseg(),
    vstrMovingSurf(),
    strMovingMri(),
    strOutput("out.mgz"),
    strOutputField("out_field.mgz"),
    strOutputMesh(),
    strOutputSurf(),        // just a placeholder
    strOutputSurfAffine(),  // just a placeholder
    strGcam(),
    strOutputAffine()  
{
  eltVolMin = 2;
  eltVolMax = 21;
  poissonRatio = .3f;
  YoungModulus = 10;
  iSteps = 1;    // by default, simple linear elastic model
  iEndStep = -1; // by default, do an extra step to finish converging
  surfSubsample = -1;
  bUseOldTopologySolver = false;
  bUsePialForSurf = false;
}

int
IoParams::parse(int argc, char* argv[], std::string& errMsg)
{
  // Define long options for getopt_long
  static struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"fixed-mri", required_argument, 0, 0},
    {"moving-mri", required_argument, 0, 0},
    {"aseg", required_argument, 0, 0},
    {"fixed-surf", required_argument, 0, 0},
    {"moving-surf", required_argument, 0, 0},
    {"aparc", required_argument, 0, 0},
    {"out", required_argument, 0, 'o'},
    {"out-field", required_argument, 0, 0},
    {"out-mesh", required_argument, 0, 0},
    {"out-surf", required_argument, 0, 0},
    {"out-surf-affine", required_argument, 0, 0},
    {"gcam", required_argument, 0, 0},
    {"out-affine", required_argument, 0, 0},
    {"dbg-output", required_argument, 0, 0},
    {"elt-vol", required_argument, 0, 0},
    {"elt-vol-min", required_argument, 0, 0},
    {"elt-vol-max", required_argument, 0, 0},
    {"poisson", required_argument, 0, 0},
    {"young", required_argument, 0, 0},
    {"surf-subsample", required_argument, 0, 0},
    {"cache-transform", required_argument, 0, 0},
    {"compress", no_argument, 0, 0},
    {"fem-steps", required_argument, 0, 0},
    {"fem-end-step", required_argument, 0, 0},
    {"topology-old", no_argument, 0, 0},
    {"use-pial-for-surf", no_argument, 0, 0},
    {"penalty-weight", required_argument, 0, 0},
    {0, 0, 0, 0}
  };

  int option_index = 0;
  int c;

  // Counter for multiple surfaces
  int fixed_surf_count = 0;
  int moving_surf_count = 0;
  int aparc_count = 0;

  while ((c = getopt_long(argc, argv, "ho:", long_options, &option_index)) != -1)
  {
    if (c == 0)
    {
      // Long option
      std::string opt_name = long_options[option_index].name;

      if (opt_name == "fixed-mri")
        strFixedMri = optarg;
      else if (opt_name == "moving-mri")
        strMovingMri = optarg;
      else if (opt_name == "aseg")
        strAseg = optarg;
      else if (opt_name == "fixed-surf")
        vstrFixedSurf.push_back(optarg);
      else if (opt_name == "moving-surf")
        vstrMovingSurf.push_back(optarg);
      else if (opt_name == "aparc")
      {
        hasAparc = true;
        vstrAparc.push_back(optarg);
      }
      else if (opt_name == "out-field")
        strOutputField = optarg;
      else if (opt_name == "out-mesh")
        strOutputMesh = optarg;
      else if (opt_name == "out-surf")
        strOutputSurf = optarg;
      else if (opt_name == "out-surf-affine")
        strOutputSurfAffine = optarg;
      else if (opt_name == "gcam")
        strGcam = optarg;
      else if (opt_name == "out-affine")
        strOutputAffine = optarg;
      else if (opt_name == "dbg-output")
        strDebug = optarg;
      else if (opt_name == "elt-vol")
        eltVolMin = eltVolMax = atof(optarg);
      else if (opt_name == "elt-vol-min")
        eltVolMin = atof(optarg);
      else if (opt_name == "elt-vol-max")
        eltVolMax = atof(optarg);
      else if (opt_name == "poisson")
        poissonRatio = atof(optarg);
      else if (opt_name == "young")
        YoungModulus = atof(optarg);
      else if (opt_name == "surf-subsample")
        surfSubsample = atof(optarg);
      else if (opt_name == "cache-transform")
        strTransform = optarg;
      else if (opt_name == "compress")
        compress = true;
      else if (opt_name == "fem-steps")
        iSteps = atoi(optarg);
      else if (opt_name == "fem-end-step")
        iEndStep = atoi(optarg);
      else if (opt_name == "topology-old")
        bUseOldTopologySolver = true;
      else if (opt_name == "use-pial-for-surf")
        bUsePialForSurf = true;
      else if (opt_name == "penalty-weight")
        ; // This will be handled by solver directly, but accept it here
    }
    else if (c == 'h')
    {
      help_exit();
    }
    else if (c == 'o')
    {
      strOutput = optarg;
    }
    else if (c == '?')
    {
      // getopt_long already printed an error message
      errMsg += "Invalid command-line option.\n";
      return 1;
    }
  }

  // Validate required parameters
  if (strFixedMri.empty())
    errMsg += " No fixed volume specified (use --fixed-mri)\n";

  if (strMovingMri.empty())
    errMsg += " No moving volume specified (use --moving-mri)\n";

  if (vstrFixedSurf.empty())
    errMsg += " No fixed surface specified (use --fixed-surf)\n";

  if (vstrMovingSurf.empty())
    errMsg += " No moving surface specified (use --moving-surf)\n";

  // Check surface count match
  if (vstrFixedSurf.size() != vstrMovingSurf.size())
  {
    errMsg += " Number of fixed and moving surfaces must match\n";
  }

  // Check aparc count if provided
  if (hasAparc && vstrAparc.size() != vstrFixedSurf.size())
  {
    errMsg += " Number of aparc files must match number of surfaces\n";
  }

  return 0;
}

/*
Original PETSc-based parse function has been removed during MFEM migration.
Below is the outline of what it parsed for reference:

Required: -fixed_mri, -moving_mri, -fixed_surf, -moving_surf
Optional: -aseg, -aparc, -out, -out_field, -out_mesh, -out_surf,
          -out_surf_affine, -gcam, -out_affine, -dbg_output,
          -elt_vol, -elt_vol_range, -poisson, -young, -surf_subsample,
          -cache_transform, -compress, -fem_steps, -fem_end_step,
          -topology_old, -use_pial_for_surf

TODO: Implement new command-line parser using getopt or similar library.
*/


void
IoParams::help_exit()
{
  std::cout << "\nUsage: surf2vol [OPTIONS]\n\n"
  << "Diffuse surface deformation to volumes using elastic finite element methods.\n"
  << "Migrated from PETSc to MFEM library.\n\n"
  << "Required arguments:\n"
  << "  --fixed-mri <file>        Fixed volume MRI file\n"
  << "  --moving-mri <file>       Moving volume MRI file\n"
  << "  --fixed-surf <file>       Fixed surface file (can be specified multiple times)\n"
  << "  --moving-surf <file>      Moving surface file (can be specified multiple times)\n"
  << "\nOptional arguments:\n"
  << "  -h, --help                Show this help message\n"
  << "  --aseg <file>             ASEG file for fixed volume\n"
  << "  --aparc <file>            APARC file (can be specified multiple times, must match surfaces)\n"
  << "  -o, --out <file>          Output file (default: out.mgz)\n"
  << "  --out-field <file>        Output field file (default: out_field.mgz)\n"
  << "  --out-mesh <file>         Output mesh file (.tm3d)\n"
  << "  --out-surf <file>         Output surface file root\n"
  << "  --out-surf-affine <file>  Output affine surface file root\n"
  << "  --gcam <file>             Output GCAM file\n"
  << "  --out-affine <file>       Output affine morphed volume\n"
  << "  --dbg-output <file>       Debug output (writes morph at each iteration)\n"
  << "\nFEM parameters:\n"
  << "  --elt-vol <value>         Element volume (sets both min and max)\n"
  << "  --elt-vol-min <value>     Minimum element volume (default: 2)\n"
  << "  --elt-vol-max <value>     Maximum element volume (default: 21)\n"
  << "  --poisson <value>         Poisson ratio (default: 0.3, must be < 0.5)\n"
  << "  --young <value>           Young's modulus (default: 10)\n"
  << "  --fem-steps <n>           Number of FEM steps (default: 1)\n"
  << "  --fem-end-step <n>        End step (default: -1)\n"
  << "  --penalty-weight <value>  Penalty weight for MFC conditions\n"
  << "\nOther options:\n"
  << "  --surf-subsample <value>  Surface subsampling distance\n"
  << "  --cache-transform <file>  Cache/load transform from file\n"
  << "  --compress                Compress morph at each step\n"
  << "  --topology-old            Use old topology solver\n"
  << "  --use-pial-for-surf       Use pial surface for mesh construction\n"
  << "\nExample:\n"
  << "  surf2vol --fixed-mri fixed.mgz --moving-mri moving.mgz \\\n"
  << "           --fixed-surf lh.white --moving-surf lh.white.moved \\\n"
  << "           --out output.mgz\n\n";
  exit(1);
}

/*-----------------------------------

File input

-----------------------------------*/

class SurfaceMarker
{
public:
  SurfaceMarker(MRIS* psurf)
  {
    m_psurf = psurf;

    // BECAUSE OF THE VERY REASON THE LABELING IS SOMEWHAT ARBITRARY,
    //
    // IT DOES NOT MEAN THE SURFACE REGISTRATION IS BAD IN THAT AREA
    // HENCE DO NOT DISCARD INFORMATION IN BANKSSTS OR FRONTALPOLE
    this->addLabelToIgnored("unknown");
    //this->addLabelToIgnored("bankssts");
    this->addLabelToIgnored("corpuscallosum");
    //this->addLabelToIgnored("frontalpole");

    std::sort(m_vecAparcIgnored.begin(),
              m_vecAparcIgnored.end() );

    std::cout << " ignored labels = ";
    std::copy(m_vecAparcIgnored.begin(),
              m_vecAparcIgnored.end(),
              std::ostream_iterator<int>(std::cout, " " ) );
    std::cout << "\n";
  }
  void operator()(VERTEX* pvtx)
  {
    int ctIndex;
    CTABfindAnnotation( m_psurf->ct, pvtx->annotation, &ctIndex );
    if ( std::binary_search( m_vecAparcIgnored.begin(),
                             m_vecAparcIgnored.end(),
                             ctIndex )
       )
      pvtx->ripflag = 1;
  }
private:
  SurfaceMarker()
{};

  std::vector<int> m_vecAparcIgnored;
  MRIS* m_psurf;

  void
  addLabelToIgnored(const char* name)
  {
    int aparcVal;
    CTABfindName( m_psurf->ct, const_cast<char*>(name), &aparcVal );
    std::cout << " label = " << name << " -> value = " << aparcVal << std::endl;
    m_vecAparcIgnored.push_back(aparcVal);
  }
};

static
int
process_input_data( const IoParams& params,
                    MRI* &mri_fixed, SurfaceVectorType& vmris_fixed,
                    MRI* &mri_moving, SurfaceVectorType& vmris_moving,
                    MRI* &mri_aseg )
{
  MRI_SURFACE* mris;

  // process input volumes
  std::cout << " reading fixed volume\n";
  mri_fixed = MRIread( const_cast<char*>( params.strFixedMri.c_str()) );
  if ( !mri_fixed )
  {
    std::cerr << " Error reading volume file " << params.strFixedMri
    << std::endl;
    exit(1);
  }

  std::cout << " reading moving volume\n";
  mri_moving = MRIread( const_cast<char*>( params.strMovingMri.c_str()) );
  if ( !mri_moving )
  {
    std::cerr << " Error reading volume file " << params.strMovingMri
    << std::endl;
    exit(1);
  }

  if ( !params.strAseg.empty() )
  {
    std::cout << " reading fixed MRI aseg labels\n";
    mri_aseg = MRIread( const_cast<char*>(params.strAseg.c_str()) );
    if ( !mri_aseg )
    {
      std::cerr << " Error reading aseg volume file "
      << params.strAseg << std::endl;
      exit(1);
    }
  }

  for ( IoParams::StringVectorType::const_iterator
        cit = params.vstrFixedSurf.begin();
        cit != params.vstrFixedSurf.end(); ++cit )
  {
    mris = MRISread( const_cast<char*>( cit->c_str() ) );
    if (!mris)
    {
      std::cerr << " Error reading surface "
      << *cit << std::endl;
      exit(1);
    }
    SurfacePointer sp;
    sp.GetTypeFromName(*cit);
    sp.mris = mris;

    vmris_fixed.push_back( sp );
  } // next cit

  if ( params.hasAparc )
  {
    std::cout << " processing APARC\n";
    // the aparc vector must have the same size as the fixed surf vector
    if ( params.vstrAparc.size() != vmris_fixed.size() )
      std::cerr << " Size mismatch between APARC vector and vmris_fixed\n";
    else
    {
      IoParams::StringVectorType::const_iterator citAparc
      = params.vstrAparc.begin();
      for ( SurfaceVectorType::iterator it = vmris_fixed.begin();
            it != vmris_fixed.end();
            ++it, ++citAparc )
      {
        // allow skipping one aparc
        if ( *citAparc == "none" ) continue;

        if ( MRISreadAnnotation( it->mris,
                                 const_cast<char*>(citAparc->c_str()) )
             != NO_ERROR )
        {
          std::cerr << " Error reading APARC file " << *citAparc
          << std::endl;
          exit(1);
        }

        SurfaceVertexIterator  vtxIter(it->mris);
        SurfaceMarker functor( it->mris );
        vtxIter.Execute( functor );

      } // next it, citAparc
    }
  }

  for ( IoParams::StringVectorType::const_iterator
        cit = params.vstrMovingSurf.begin();
        cit != params.vstrMovingSurf.end(); ++cit )
  {
    mris = MRISread( const_cast<char*>( cit->c_str() ) );
    if (!mris)
    {
      std::cerr << " Error reading surface "
      << *cit << std::endl;
      exit(1);
    }
    SurfacePointer sp;
    sp.GetTypeFromName(*cit);
    sp.mris = mris;
    vmris_moving.push_back(sp);
  } // next cit

  return 0;
}

//===================================================

void output_bc_locations( const tSolver& solver,
                          MRI* mri)
{
  typedef tSolver::BcContainerConstIterator BcConstIterator;
  BcConstIterator begin, end;

  solver.getBcIterators(begin,end);

  MRI* mri_out = MRIalloc( mri->width,
                           mri->height,
                           mri->depth,
                           MRI_FLOAT);
  MRIcopyHeader( mri, mri_out);
  MRIvalueFill( mri_out, 0.0f);

  int x,y,z;
  for ( ; begin!=end; ++begin)
  {
    x = (int)std::floor( (*begin)->pt(0) + .5 );
    y = (int)std::floor( (*begin)->pt(1) + .5 );
    z = (int)std::floor( (*begin)->pt(2) + .5 );

    MRIsetVoxVal( mri_out, x,y,z,0, 250.0f);
  }

  MRIwrite( mri_out, const_cast<char*>("bc.mgz") );
}

void
compute_bc_error( const tSolver& solver )
{
  double dActiveError = 0.0;
  double dTotalError = 0.0;
  double dcount = 0;

  typedef tSolver::BcContainerConstIterator BcConstIterator;
  BcConstIterator it, end;

  solver.getBcIterators(it,end);



  for ( ; it != end; ++it )
  {
    double dNorm = ( solver.get_mesh()->dir_img( (*it)->pt )
                     - (*it)->pt - (*it)->delta ).norm();
    if ( (*it)->isActive )
      dActiveError += dNorm;
    dTotalError += dNorm;
    dcount += 1.0;
  } // next it

  if ( dcount == 0 )
  {
    std::cerr << " NO BCs !?!?\n";
    return;
  }

  std::cout << " Active Residual error = " << dActiveError/dcount << std::endl
  << " Total error = " << dTotalError / dcount << std::endl;

}


#ifdef USE_SURF_FORWARD_MORPH
static void surf_forward_morph(const CMesh3d* pmesh, std::string strName,
                               SurfaceVectorType& vmris_moving,
                               MRI* mri_fixed, TCoords<double,3> cmin,
                               TCoords<double,3> cmax)
{

  CMesh3d meshInverse( *(&(*pmesh)));
  ;

  meshInverse.invert();
  meshInverse.build_index_src();

  // create the transforms
  gmp::FemTransform3d femTransform;
  femTransform.set_mesh( std::shared_ptr<TMesh3d>(&meshInverse) );

  femTransform.m_signalTopology = true;

  for ( SurfaceVectorType::iterator it = vmris_moving.begin();
        it != vmris_moving.end(); ++it )
  {
    std::cout << " processing surface "
    << it - vmris_moving.begin()
    << std::endl;
    MRI_SURFACE* mris = it->mris;
    VERTEX* pvtx = &( mris->vertices[0] );
    unsigned int nvertices = (unsigned int)mris->nvertices;

    int errCount = 0;
    std::cout << " bbox for inverse stuff = " << cmin << " , " << cmax << std::endl;
    gmp::FemTransform3d::tCoords img;
    for (unsigned int ui=0;
         ui < nvertices;
         ++ui, ++pvtx )
    {
      gmp::FemTransform3d::tCoords pt;
      pt(0) = pvtx->x;
      pt(1) = pvtx->y;
      pt(2) = pvtx->z;

      img = femTransform.img( pt );

      if ( !img.isValid() )
      {
        std::cerr << " Invalid FemTransform at position " << pt << std::endl;
        ++errCount;
        if ( errCount > 100 ) exit(1);
        continue;
      }
      if ( img(0)<0 || img(0)>mri_fixed->width ||
           img(1)<0 || img(1)>mri_fixed->height ||
           img(2)<0 || img(2)>mri_fixed->depth )
        continue;

      pvtx->x = img(0);
      pvtx->y = img(1);
      pvtx->z = img(2);

    } // next ui , pvtx

    // save the surface
    char fname[256];
    sprintf(fname,"%s_%d.white", strName.c_str(),
            (int)(it-vmris_moving.begin()) );
    std::cout << " writing surface to file " << fname
    << std::endl;
    convert_vox_to_surf(mris, mri_fixed);
    MRISwrite( mris, fname);
  } // next it

}
#endif //ifdef USE_SURF_FORWARD_MORPH


static void compute_fem_error(const Transform3SPointer ptransform,
                              const PointsContainerType& container)
{
  std::cout << " FEM error\n";

  double dInitial(.0), dFinal(.0), count(.0);
  double dMaxError(.0), dCrtError;

  int invalid = 0;
  for ( PointsContainerType::const_iterator cit = container.begin();
        cit != container.end();
        ++cit)
  {
    // get the image of the point
    TCoords<double,3> img = ptransform->img( cit->first );
    if ( !img.isValid() )
    {
      ++invalid;
      continue;
    }

    dCrtError = ( img - cit->second ).norm();
    dFinal += dCrtError;
    dInitial += (cit->second - cit->first).norm();
    count += 1.0;
    dMaxError = std::max( dMaxError, dCrtError );
  } // next cit

  if (!count)
  {
    std::cerr << " Error with count\n";
    return;
  }

  std::cout << " Initial displacement = " << dInitial/count << std::endl
  << " Final displacement = " << dFinal/count << std::endl
  << " Max after displacement = " << dMaxError << std::endl
  << " Error points = " << invalid << std::endl;
}

#ifdef USE_GATHER_DEFECT_INFO
static void gather_defect_info( CMesh3d::ElementIndexContainer& lstElts,
                                tSolver::BcMfcInfoType& mfcInfo)
{
  int bcDefects = 0;
  tSolver::BcMfcInfoType::const_iterator bcIter;
  for ( CMesh3d::ElementIndexContainer::const_iterator cit = lstElts.begin();
        cit != lstElts.end(); ++cit )
  {
    bcIter = mfcInfo.find( *cit );
    if ( bcIter != mfcInfo.end() )
    {
      bcDefects++;
      std::cout << " elt index " << *cit << " -> norm var = "
      << bcIter->second.second << std::endl;
    }
  } // next cit

  std::cout << " bc defects = " << bcDefects << std::endl;
}
#endif //ifdef USE_GATHER_DEFECT_INFO

tDblCoords cross_product(const tDblCoords& a,
                         const tDblCoords& b)
{
  tDblCoords result;

  result(0) = a(1)*b(2) - a(2)*b(1);
  result(1) = a(2)*b(0) - a(0)*b(2);
  result(2) = a(0)*b(1) - a(1)*b(0);

  return result;
}

void check_surface_defects(const SurfaceVectorType& vmris_x,
                           const SurfaceVectorType& vmris_fx)
{
  SurfaceVectorType::const_iterator cit_x, cit_fx;

  for ( cit_x = vmris_x.begin(), cit_fx = vmris_fx.begin();
        cit_x != vmris_x.end();
        ++cit_x, ++cit_fx)
  {
    MRI_SURFACE* mrisX = cit_x->mris;
    MRI_SURFACE* mrisFx = cit_fx->mris;

    MRIScomputeMetricProperties( mrisX );
    MRIScomputeMetricProperties( mrisFx );

    FACE* faceX = &( mrisX->faces[0] );
    FACE* faceFx = &( mrisFx->faces[0] );
    unsigned int nfaces = (unsigned int)( mrisX->nfaces );

    for (unsigned int ui = 0;
         ui < nfaces;
         ++ui, ++faceX, ++faceFx )
    {
      if ( faceX->area * faceFx->area < 0.0f )
        std::cout << " found defect for face " << ui << std::endl;
    }
  } // next cit_x, cit_fx
}

int nint(float f)
{
  return (int)std::floor(f+.5f);
}

void
dbg_surf2vol(const PointsContainerType& c,
             MRI* mri)
{
  MRI* out_fixed = MRIalloc( mri->width,
                             mri->height,
                             mri->depth,
                             MRI_UCHAR);
  MRIcopyHeader( mri, out_fixed );
  MRIvalueFill(out_fixed, 0);

  MRI* out_moving = MRIalloc( mri->width,
                              mri->height,
                              mri->depth,
                              MRI_UCHAR);
  MRIcopyHeader( mri, out_moving);
  MRIvalueFill(out_moving, 0);

  for ( PointsContainerType::const_iterator cit = c.begin();
        cit != c.end(); ++cit )
  {
    MRIvox(out_fixed,
           nint( cit->first(0) ),
           nint( cit->first(1) ),
           nint( cit->first(2) ) ) = 200;
    MRIvox(out_moving,
           nint( cit->second(0) ),
           nint( cit->second(1) ),
           nint( cit->second(2) ) ) = 200;
  } // next cit

  MRIwrite( out_fixed, "dbg_fixed.mgz");
  MRIwrite( out_moving, "dbg_moving.mgz");

  MRIfree(&out_fixed);
  MRIfree(&out_moving);

}

unsigned int
update_bc_container(PointsContainerType& container,
                    Transform3SPointer ptransform)
{
  unsigned int countErased(0);
  for ( PointsContainerType::iterator
        it = container.begin();
        it != container.end();
      )
  {
    it->first = ptransform->img(  it->first );
    if ( !it->first.isValid() )
    {
      it = container.erase(it);
      ++countErased;
    }
    else ++it;
  }
  return countErased;
}
