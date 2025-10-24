
#ifndef H_SOLVER_H
#define H_SOLVER_H

#include <assert.h>

#include <functional>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>

#include "mfem.hpp"

#include "mesh.h"
#include "cstats.h"

//----------------------------------------------------
//
// class declaration
//
//----------------------------------------------------

template<class Cstr, int n>
struct BC
{
  typedef TCoords<double,n> tCoords;
  typedef TMesh<Cstr,n> tMesh;

  tCoords pt;
  tCoords delta;
  bool    isActive;
  BC() : pt(), delta(), isActive(false)
  {}
  BC(const tCoords& _pt, const tCoords& _delta) : pt(_pt), delta(_delta),
      isActive(false)
  {}
  virtual ~BC()
  {};
  virtual bool find_candidate(tMesh* pmesh)
  {
    return false;
  }
};

template<class Cstr, int n>
struct BCNatural : public BC<Cstr,n>
{
  typedef typename BC<Cstr, n>::tCoords tCoords;
  typedef typename BC<Cstr,n>::tMesh tMesh;
  typedef TNode<n> tNode;

  tNode* pnode;
  BCNatural() : BC<Cstr, n>(), pnode(NULL)
  {}
  BCNatural(const tCoords& _pt, const tCoords& _delta) : BC<Cstr,n>(_pt, _delta), pnode(NULL)
  {}
  virtual ~BCNatural()
  {};
  virtual bool find_candidate(tMesh* pmesh)
  {
    pnode = pmesh->closest_node(this->pt);
    return (pnode!=NULL);
  }
};

template<class Cstr, int n>
struct BCMfc : public BC<Cstr, n>
{
  typedef typename BC<Cstr,n>::tCoords tCoords;
  typedef TElement<n> tElement;
  typedef typename BC<Cstr,n>::tMesh tMesh;

  tElement* pelt;
  BCMfc(const tCoords& _pt, const tCoords& _delta) : BC<Cstr,n>(_pt, _delta), pelt(NULL)
  {}
  virtual ~BCMfc()
  {};
  virtual bool find_candidate(tMesh* pmesh)
  {
    pelt = pmesh->element_at_point(this->pt);
    return (pelt!=NULL);
  }
};

//--------------------------------------------------------

template<class Cstr,int n>
class TSolver
{
public:
  typedef TCoords<int,n> tIntCoords;
  typedef TCoords<double,n> tCoords;
  typedef TMesh<Cstr,n> tMesh;
  typedef TNode<n> tNode;
  typedef TElement<n> tElement;

  typedef BC<Cstr,n> tBC;
  typedef BCNatural<Cstr,n> tBCNatural;
  typedef BCMfc<Cstr,n> tBCMfc;

  TSolver();
  TSolver(tIntCoords&);
  virtual ~TSolver();

  void clone(const TSolver&);
  void add_bc_natural(const tCoords& pt,
                      const tCoords& delta);
  void add_bc_natural(tNode* pnode,
                      const  tCoords& delta);
  void add_bc_mfc(const tCoords& pt,
                  const tCoords& delta);

  virtual int solve();
  void set_mesh(tMesh* pmesh)
  {
    m_pmesh = pmesh;
  }
  const tMesh* get_mesh() const
  {
    return m_pmesh;
  }

  void set_displayLevel(int level)
  {
    m_displayLevel = level;
  }

  int get_displayLevel() const
  {
    return m_displayLevel;
  }

  void setThreshold(double dval)
  {
    m_bcThreshold = dval;
    m_useThreshold=true;
  }

  // Advanced solver configuration methods
  void setSolverTolerance(double tol) { m_solverTolerance = tol; }
  void setSolverMaxIter(int maxIter) { m_solverMaxIter = maxIter; }
  void setSolverType(const std::string& type) { m_solverType = type; }
  void setPreconditioner(const std::string& pc) { m_precondType = pc; }
  void setInitialGuessNonzero(bool flag) { m_initGuessNonzero = flag; }
  void setDebugPrint(const std::string& path) { m_debugPrintPath = path; }

  typedef std::vector<tBC*> BcContainerType;
  typedef typename BcContainerType::const_iterator BcContainerConstIterator;
  unsigned int getBcIterators(BcContainerConstIterator& begin,
                              BcContainerConstIterator& end) const
  {
    begin = m_vBc.begin();
    end = m_vBc.end();
    return m_vBc.size();
  }

  double m_mfcWeight;

  int check_bc_error(double& dRemainingRatio);

  // mainly for debugging purposes
  // when assigning BC MFC - keep the information about
  // the element available for later probing
  typedef typename
  std::map<unsigned int, std::pair<tCoords, double> > BcMfcInfoType;
  BcMfcInfoType m_mfcInfo;


protected:
  BcContainerType m_vBc;

  mfem::SparseMatrix*   m_stiffness;
  mfem::Vector*   m_load;
  mfem::Vector*   m_delta;

  tMesh* m_pmesh;

  int  done_bc_natural(); // distribute the BC
  int  done_bc_mfc();
  int  setup_matrix(bool showInfo=false); // assembly the stiffness matrix
  int  setup_load(); // assembly force load by reduction of the LHS
  int  setup_load_sym(); // assembly force load by
  // conditioning both rows and cols
  int  setup_load_mfc();
  int  comm_solution(); // set sol values in respective nodes

  int  add_elt_matrix(const tElement* pelt);

  int  add_elt_mfc_lhs(const tElement* pelt, tCoords& pt);
  int  add_elt_mfc_rhs(const tElement* pelt, tCoords& pt, tCoords& delta);

  int  m_displayLevel; // 0=critical, 1=important, 2=detailed

  bool m_useThreshold; // sets whether a threshold should
  // be used or not when setting the BC
  double m_bcThreshold;

  // Advanced solver configuration
  double m_solverTolerance;      // Convergence tolerance
  int    m_solverMaxIter;        // Maximum iterations
  std::string m_solverType;      // Solver type: "cg", "gmres", "minres"
  std::string m_precondType;     // Preconditioner: "none", "jacobi", "gs"
  bool   m_initGuessNonzero;     // Use nonzero initial guess
  std::string m_debugPrintPath;  // Path for debug matrix output
};


//--------------------------------
//
// Derived class
// uses a direct solver instead of KSP
//
//  only handles natural BCs
//---------------------------------

template<class Cstr, int n>
class TDirectSolver : public TSolver<Cstr, n>
{
public:
  TDirectSolver();
  ~TDirectSolver();

  int solve();
};

//--------------------------------------------------------------------
//
// class implementation
//
//--------------------------------------------------------------------

template<class Cstr,int n>
TSolver<Cstr,n>::TSolver()
{
  m_pmesh = NULL;

  m_stiffness = 0;
  m_load = 0;
  m_delta = 0;

  m_displayLevel = 1;

  m_useThreshold = false;
  m_bcThreshold = 0.0;

  m_mfcWeight = 1.0;

  // Initialize advanced solver options with defaults
  m_solverTolerance = 1.0e-9;
  m_solverMaxIter = 10000;
  m_solverType = "cg";        // Default to Conjugate Gradient
  m_precondType = "none";     // No preconditioner by default
  m_initGuessNonzero = false;
  m_debugPrintPath = "";      // No debug output by default
}

template<class Cstr, int n>
TSolver<Cstr,n>::~TSolver()
{
  for ( typename BcContainerType::iterator it = m_vBc.begin();
        it != m_vBc.end(); ++it )
    delete *it;
  m_vBc.clear();
}

template<class Cstr,int n>
void
TSolver<Cstr,n>::clone(const TSolver& s)
{
  m_vBc   = s.m_vBc;
  //m_ticks = s.m_ticks;
  m_pmesh  = s.m_pmesh;

  m_stiffness = 0;
  m_load = 0;
  m_delta = 0;
}

template<class Cstr,int n>
void
TSolver<Cstr,n>::add_bc_natural(const tCoords& pt,
                                const tCoords& delta)
{
  tBCNatural *bc = new tBCNatural(pt,delta);
  m_vBc.push_back(bc);
}

template<class Cstr, int n>
void
TSolver<Cstr,n>::add_bc_natural(tNode* pnode,
                                const tCoords& delta)
{
  tBCNatural *bc = new tBCNatural(pnode->coords(),
                                  delta );
  bc->pnode = pnode;
  bc->pt = pnode->coords();
  m_vBc.push_back(bc);
}

template<class Cstr, int n>
void
TSolver<Cstr,n>::add_bc_mfc(const tCoords& pt,
                            const tCoords& delta)
{
  tBCMfc *bc = new tBCMfc(pt, delta);
  m_vBc.push_back(bc);
}

#undef __FUNCT__
#define __FUNCT__ "TSolver::solve"
template<class Cstr,int n>
int
TSolver<Cstr,n>::solve()
{
  std::cout << " penalty_weight = " << m_mfcWeight << std::endl;

  done_bc_natural();
  done_bc_mfc();

  setup_matrix(true); // display information by default for the major system

  // MFEM matrices are finalized during construction - no explicit assembly needed

  setup_load_sym(); // pin down Natural conditions
  setup_load_mfc(); // setup MFC conditions

  // MFEM vectors are ready to use after setup

  // Debug: Print matrix if requested
  if (!m_debugPrintPath.empty())
  {
    std::ofstream matFile(m_debugPrintPath + "_stiffness.txt");
    m_stiffness->Print(matFile);
    matFile.close();
    std::cout << " Wrote stiffness matrix to " << m_debugPrintPath << "_stiffness.txt" << std::endl;

    std::ofstream vecFile(m_debugPrintPath + "_load.txt");
    m_load->Print(vecFile);
    vecFile.close();
    std::cout << " Wrote load vector to " << m_debugPrintPath << "_load.txt" << std::endl;
  }

  // Allocate solution vector
  m_delta = new mfem::Vector(m_load->Size());
  if (m_initGuessNonzero)
  {
    *m_delta = *m_load; // Use load as initial guess
    std::cout << " Using nonzero initial guess for system solving\n";
  }
  else
  {
    *m_delta = 0.0; // Initialize to zero
  }

  // Create preconditioner if requested
  mfem::Solver* precond = nullptr;
  if (m_precondType == "jacobi" || m_precondType == "gs")
  {
    int sweeps = 1;
    if (m_precondType == "jacobi")
    {
      precond = new mfem::DSmoother(0, 1.0, sweeps); // Jacobi (type 0)
      std::cout << " Using Jacobi preconditioner (" << sweeps << " sweep)\n";
    }
    else if (m_precondType == "gs")
    {
      precond = new mfem::DSmoother(1, 1.0, sweeps); // Gauss-Seidel (type 1)
      std::cout << " Using Gauss-Seidel preconditioner (" << sweeps << " sweeps)\n";
    }
    precond->SetOperator(*m_stiffness);
  }

  // Create and configure the iterative solver based on type
  mfem::IterativeSolver* solver = nullptr;

  if (m_solverType == "gmres")
  {
    mfem::GMRESSolver* gmres = new mfem::GMRESSolver();
    gmres->SetKDim(50); // Restart parameter
    solver = gmres;
    std::cout << " Using GMRES solver (restart=" << 50 << ")\n";
  }
  else if (m_solverType == "minres")
  {
    solver = new mfem::MINRESSolver();
    std::cout << " Using MINRES solver\n";
  }
  else // default to CG
  {
    solver = new mfem::CGSolver();
    std::cout << " Using Conjugate Gradient solver\n";
  }

  // Configure solver parameters
  solver->SetRelTol(m_solverTolerance);
  solver->SetMaxIter(m_solverMaxIter);
  solver->SetPrintLevel(m_displayLevel > 1 ? 2 : 0);
  solver->SetOperator(*m_stiffness);

  if (precond)
  {
    solver->SetPreconditioner(*precond);
  }

  std::cout << " Solver configuration: tol=" << m_solverTolerance
            << ", maxiter=" << m_solverMaxIter << std::endl;

  // Solve the linear system A*x = b
  solver->Mult(*m_load, *m_delta);

  if (solver->GetConverged())
  {
    std::cout << " MFEM " << m_solverType << " solver converged in "
              << solver->GetNumIterations()
              << " iterations with final norm " << solver->GetFinalNorm() << std::endl;
  }
  else
  {
    std::cout << " WARNING: MFEM " << m_solverType
              << " solver did not converge after " << solver->GetNumIterations()
              << " iterations!" << std::endl;
  }

  // Check the error: ||A*x - b||
  mfem::Vector residual(m_load->Size());
  m_stiffness->Mult(*m_delta, residual);
  residual -= *m_load;
  double norm = residual.Norml2();
  std::cout << "Absolute-Norm of error = " << norm << std::endl;

  // Debug: Print solution if requested
  if (!m_debugPrintPath.empty())
  {
    std::ofstream solFile(m_debugPrintPath + "_solution.txt");
    m_delta->Print(solFile);
    solFile.close();
    std::cout << " Wrote solution vector to " << m_debugPrintPath << "_solution.txt" << std::endl;
  }

  // Clean up
  delete solver;
  if (precond) delete precond;

  //--------
  comm_solution(); // distribute obtained displacements to nodes

  double dremainingRatio;
  check_bc_error(dremainingRatio);

  // release MFEM resources
  delete m_delta;
  m_delta = nullptr;
  delete m_load;
  m_load = nullptr;
  delete m_stiffness;
  m_stiffness = nullptr;

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "TSolver::done_bc_natural"
template<class Cstr,int n>
int
TSolver<Cstr,n>::done_bc_natural()
{
  std::vector<tCoords> vdelta; // holds the point-wise diff for each BC
  bool bFailed = false;

  tNode* pnode = NULL;
  for ( typename BcContainerType::iterator it = m_vBc.begin();
        it != m_vBc.end();
        ++it )
  {
    if ( tBCNatural* bc = dynamic_cast<tBCNatural*>( *it ) )
    {
      // if node was not previously specified, find closest now
      if (!bc->pnode)
      {
        pnode = m_pmesh->closest_node(bc->pt);
        bc->pnode = pnode;
      }
      else
        pnode = bc->pnode;

      if ( pnode )
      {
        if ( !m_useThreshold ||
             (pnode->coords()- bc->pt).norm() < m_bcThreshold )
        {
          pnode->set_bc(bc->delta);
          vdelta.push_back( pnode->coords() - bc->pt );

          bc->isActive = true;

          if ( m_displayLevel>1 )
            std::cout << "setting bc " << pnode->coords()
            << " -> " << bc->delta << "\n"
            << "\t instead " << bc->pt << " -> norm = " << vdelta.back().norm()
            << std::endl;
        }
      }
      else
      {
        std::cerr
        << "TSolver::done_bc_natural -> failed to find node close to "
        << bc->pt << std::endl;
        bFailed = true;
      }
    }
  }

  if ( m_displayLevel &&
       bFailed ) std::cout << " !!!!! There were FAILED BCs\n";
  if ( m_displayLevel )
  {
    std::cout
    <<  " computing statistics for the displacement application error\n";
    double dAvgNorm = 0.0;
    for ( typename std::vector< tCoords>::const_iterator cit = vdelta.begin();
          cit != vdelta.end();
          ++cit )
      dAvgNorm += cit->norm();
    dAvgNorm /= (double)vdelta.size();
    std::cout
    << " average norm of error in placement = " << dAvgNorm << std::endl;
  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "TSolver::done_bc_mfc"
template<class Cstr, int n>
int
TSolver<Cstr,n>::done_bc_mfc()
{
  bool bFailed = false;

  typedef std::map< int, std::vector<int> > MfcCandidateType;
  MfcCandidateType candidates;
  MfcCandidateType::iterator mapIter;

  tElement* pelt = NULL;
  int index = 0;
  std::cout << " iterating\n";
  for ( typename BcContainerType::iterator it = m_vBc.begin();
        it != m_vBc.end(); ++it, ++index )
  {
    if ( tBCMfc* bc = dynamic_cast<tBCMfc*>(*it) )
    {
      pelt = m_pmesh->element_at_point(bc->pt);
      if ( pelt )
      {
        mapIter = candidates.find( pelt->get_id() );
        if ( mapIter == candidates.end() )
        {
          std::vector<int> vbuf;
          vbuf.push_back( index );
          candidates[ pelt->get_id() ] = vbuf;
        }
        else
          mapIter->second.push_back( index );
      }
      else
      {
        std::cerr << " TSolver::done_bc_mfc -> failed to find elt for coords "
        << bc->pt << std::endl;
        bFailed = true;
      }
    }
  } // next it
  std::cout << " done with candidates\n";

  // go through the assignment map and compute the 3D variances
  // per element
  tCoords mean;
  int active = 0;
  double dvarcova;
  for ( mapIter = candidates.begin();
        mapIter != candidates.end();
        ++mapIter )
  {
    std::vector<tCoords> vdelta;
    for ( std::vector<int>::const_iterator cit = mapIter->second.begin();
          cit != mapIter->second.end();
          ++cit )
    {
      vdelta.push_back( m_vBc[*cit]->delta );
    } // next cit
    // compute mean and variance per elt
    // return the norm of the covariance-matrix
    dvarcova = coords_statistics( vdelta, mean);
    m_mfcInfo[ mapIter->first ] = std::make_pair( mean, dvarcova );

    // get closest BC to the mean
    std::vector<int>::const_iterator citArgmin = mapIter->second.begin();
    double dMinDist(1000);
    double dCrtDist;

#if 0
    if ( dvarcova > 1.0 ) continue;
#endif

    for ( std::vector<int>::const_iterator cit = mapIter->second.begin();
          cit != mapIter->second.end();
          ++cit)
    {
      // need to write a routine to invert the covariance
      // matrix 3x3 - should be direct
      // use determinants, i guess
      dCrtDist = ( m_vBc[*cit]->delta - mean).norm();
      if ( dCrtDist < dMinDist )
      {
        dMinDist = dCrtDist;
        citArgmin = cit;
      }
    } // next cit

    // assign BC
    ++active;
    m_vBc[*citArgmin]->isActive = true;
    dynamic_cast<tBCMfc*>(m_vBc[*citArgmin])->pelt =
      m_pmesh->fetch_elt(mapIter->first);
  } // next mapIter

  std::cout << " Active BCs = " << active << std::endl
  << " Total BCs = " << m_vBc.size() << std::endl;

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "TSolver::setup_matrix"
template<class Cstr,int n>
int
TSolver<Cstr,n>::setup_matrix(bool showInfo)
{
  int n_eqs = m_pmesh->get_no_nodes() * n; // template par = dim
  if (showInfo)  std::cout << " no-eqs = " << n_eqs << std::endl;

  // Create MFEM SparseMatrix
  // We'll use dynamic assembly, so we don't need to specify the structure in advance
  m_stiffness = new mfem::SparseMatrix(n_eqs, n_eqs);

  int iold_val = -1;
  for ( size_t i= size_t(0);
        i<m_pmesh->get_no_elts();
        ++i )
  {
    if ( showInfo )
    {
      int idone = (int)floor( double(i) / m_pmesh->get_no_elts() * 100.0 );
      if ( idone%5==0 && idone!=iold_val )
      {
        iold_val=idone;
        std::cout << "\t percentage done= " << idone << std::endl;
      }
    }
    add_elt_matrix(m_pmesh->get_elt(i));
  }

  // Finalize the sparse matrix
  m_stiffness->Finalize();

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "TSolver::add_elt_matrix"
template<class Cstr, int n>
int
TSolver<Cstr,n>::add_elt_matrix(const tElement* pelt)
{
  int *id = new int[ pelt->no_nodes() ];
  tNode* pnode = NULL;

  for (int i= 0;
       i<pelt->no_nodes();
       ++i )
  {
    if (!pelt->get_node(i,&pnode))
    {
      std::cerr << "TSolver::add_elt_matrix -> err 1\n";
      exit(1);
    }
    id[i] = pnode->get_id();
  }

  SmallMatrix elt_matrix = pelt->get_matrix();

  int a,b;
  for (int i= 0; i<pelt->no_nodes(); ++i)
    for (int j= 0; j<pelt->no_nodes(); ++j)
      for (int k=0; k<n; ++k)
        for (int l=0; l<n; ++l)
        {
          a = n*i + k;
          b = n*j + l;

          assert(id[i]>=0);
          assert(id[j]>=0);

          int row = n*id[i]+k;
          int col = n*id[j]+l;
          double value = (double)elt_matrix(a,b);

          // Add value to sparse matrix
          m_stiffness->Add(row, col, value);
        }

  delete[] id;

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "TSolver::add_elt_mfc_lhs"
template<class Cstr, int n>
int
TSolver<Cstr,n>::add_elt_mfc_lhs(const tElement* pelt,
                                 tCoords& pt)
{
  int *id = new int[ pelt->no_nodes() ];
  tNode* pnode = NULL;

  for (int i=0; i<pelt->no_nodes(); ++i)
  {
    if ( !pelt->get_node(i,&pnode) )
    {
      std::cerr << " TSolver::add_elt_mfc -> err 2\n";
      exit(1);
    }
    id[i] = pnode->get_id();
  }

  SmallMatrix elt_matrix(n, n*pelt->no_nodes());

  SmallMatrix bufMatrix;
  for ( int i=0; i<pelt->no_nodes(); ++i)
  {
    bufMatrix.identity(n);
    bufMatrix *= pelt->shape_fct(i, pt);

    elt_matrix.set_block( bufMatrix,
                          0, i*n );
  }

  // obtain the matrix by LEFT multiplying with the transpose
  bufMatrix = elt_matrix.transpose() * elt_matrix;
  bufMatrix *= m_mfcWeight;

  int a,b;
  for (int i=0; i<pelt->no_nodes(); ++i)
    for (int j=0; j<pelt->no_nodes(); ++j)
      for (int k=0; k<n; ++k)
        for (int l=0; l<n; ++l)
        {
          a = n*i + k;
          b = n*j + l;

          int row = n*id[i] + k;
          int col = n*id[j] + l;
          double value = bufMatrix(a,b);

          // Add value to sparse matrix
          m_stiffness->Add(row, col, value);
        } // next i,j,k,l

  delete[] id;

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "TSolver::setup_load_mfc"
template<class Cstr, int n>
int
TSolver<Cstr,n>::setup_load_mfc()
{
  for ( typename BcContainerType::iterator it = m_vBc.begin();
        it != m_vBc.end(); ++it )
  {
    if ( !(*it)->isActive ) continue;

    if ( tBCMfc* bc = dynamic_cast<tBCMfc*>( &(*(*it)) ) )
    {
      add_elt_mfc_lhs( bc->pelt, bc->pt );
      add_elt_mfc_rhs( bc->pelt, bc->pt, bc->delta);
    }
  } // next it

  std::cout << " after setup_load_mfc\n";
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "TSolver::add_elt_mfc_rhs"
template<class Cstr, int n>
int
TSolver<Cstr, n>::add_elt_mfc_rhs(const tElement* pelt,
                                  tCoords& pt,
                                  tCoords& delta)
{
  tNode* pnode = NULL;

  SmallMatrix elt_matrix(n, n* pelt->no_nodes() );
  SmallMatrix bufMatrix;

  for (int i=0; i<pelt->no_nodes(); ++i)
  {
    bufMatrix.identity(n);
    bufMatrix *= pelt->shape_fct(i, pt);

    elt_matrix.set_block( bufMatrix,
                          0, i*n);
  } // next i

  elt_matrix.inplace_transpose();

  SmallMatrix mfc_rhs(n, 1);
  for (int i=0; i<n; ++i)
    mfc_rhs(i,0) = delta(i) * m_mfcWeight;

  bufMatrix = elt_matrix * mfc_rhs;

  for (int i=0; i<pelt->no_nodes(); ++i)
  {
    if ( !pelt->get_node(i, &pnode) )
    {
      std::cerr << " TSolver::add_elt_mfc_rhs -> err 3\n";
      exit(1);
    }
    for (int j=0; j<n; ++j)
    {
      int idx = n* pnode->get_id() + j;
      double value = bufMatrix( n*i + j, 0);
      (*m_load)(idx) += value; // Add to load vector
    }
  } // next i

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "TSolver::comm_solution"
template<class Cstr,int n>
int
TSolver<Cstr,n>::comm_solution()
{
  tNode* pnode = NULL;
  for (size_t i=size_t(0); i<m_pmesh->get_no_nodes(); ++i)
  {
    pnode = NULL;
    m_pmesh->get_node(i,&pnode);
    if ( !pnode )
    {
      std::cerr << "TSolver::commm_solution -> err\n";
      exit(1);
    }
    for (int j=0; j<n; ++j)
    {
      int idx = pnode->get_id()*n + j;
      pnode->set_dof_val(j, (*m_delta)(idx));
    }
  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "TSolver::setup_load"
template<class Cstr,int n>
int
TSolver<Cstr,n>::setup_load()
{
  tNode* pnode;
  int no_eqs = n * m_pmesh->get_no_nodes();

  std::map<int, double> mrhs;

  for (size_t i=size_t(0); i<m_pmesh->get_no_nodes(); ++i)
  {
    if ( !m_pmesh->get_node(i, &pnode) )
    {
      std::cerr << "TSolver::setup_load -> requested node out of range\n";
      exit(1);
    }

    for (int j=0; j<n; ++j)
    {
      if (pnode->is_dof_active(j))
        mrhs[ n* pnode->get_id() + j ] = pnode->get_dof(j);
    }
  }

  // create the RHS vector
  m_load = new mfem::Vector(no_eqs);
  *m_load = 0.0; // Initialize to zero

  // Set BC values in load vector
  for ( typename std::map<int,double>::const_iterator
        cit = mrhs.begin();
        cit!= mrhs.end();
        ++cit)
  {
    (*m_load)(cit->first) = cit->second;
  }

  // Zero out rows in stiffness matrix corresponding to BCs and set diagonal to 1
  for ( typename std::map<int,double>::const_iterator
        cit = mrhs.begin();
        cit!= mrhs.end();
        ++cit)
  {
    int row = cit->first;
    // Zero out the row
    for (int col = 0; col < no_eqs; ++col)
    {
      m_stiffness->Set(row, col, 0.0);
    }
    // Set diagonal to 1
    m_stiffness->Set(row, row, 1.0);
  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "TSolver::setup_load_sym"
template <class Cstr, int n>
int
TSolver<Cstr,n>::setup_load_sym()
{
  int no_eqs = n * m_pmesh->get_no_nodes();

  std::map<int,double> mrhs;

  if ( m_displayLevel>1) std::cout << " done_bc_natural size of container = "
    << m_vBc.size() << std::endl;

  for (typename BcContainerType::iterator it = m_vBc.begin();
       it != m_vBc.end(); ++it )
  {
    if ( !(*it)->isActive ) continue;

    if ( tBCNatural* bc = dynamic_cast<tBCNatural*>(*it) )
    {
      for (int j=0; j<n; ++j)
        mrhs[ n* bc->pnode->get_id() + j ] = bc->pnode->get_dof(j);
    }
  }
  if ( m_displayLevel>1) std::cout << " after building the map\n";

  std::cout << " LOAD size = " << mrhs.size() << std::endl;

  // create RHS vector
  m_load = new mfem::Vector(no_eqs);
  *m_load = 0.0;

  // compute RHS vector using: b' = b - A*bc_values
  mfem::Vector vecBcs(no_eqs);
  vecBcs = 0.0;

  for ( typename std::map<int,double>::const_iterator
        cit = mrhs.begin();
        cit!= mrhs.end();
        ++cit)
  {
    vecBcs(cit->first) = -cit->second;
  }

  // m_load = A * vecBcs
  m_stiffness->Mult(vecBcs, *m_load);

  // Set BC values in load vector
  for ( typename std::map<int,double>::const_iterator
        cit = mrhs.begin();
        cit!= mrhs.end();
        ++cit)
  {
    (*m_load)(cit->first) = cit->second;
  }

  // Condition the matrix symmetrically
  // Zero rows and columns for BC degrees of freedom
  for ( typename std::map<int,double>::const_iterator
        cit = mrhs.begin();
        cit!= mrhs.end();
        ++cit)
  {
    int idx = cit->first;
    // Zero out row
    for (int col = 0; col < no_eqs; ++col)
    {
      m_stiffness->Set(idx, col, 0.0);
    }
    // Zero out column
    for (int row = 0; row < no_eqs; ++row)
    {
      m_stiffness->Set(row, idx, 0.0);
    }
    // Set diagonal to 1
    m_stiffness->Set(idx, idx, 1.0);
  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "TSolver::check_bc_error"
template <class Cstr, int n>
int
TSolver<Cstr,n>::check_bc_error(double& dRemainingRatio)
{
  // go through the BCs and measure the error
  double dSum(.0), dSumConditional(.0), dSumInitial(.0), dval;
  int count(0), countConditional(0);

  tCoords img;
  unsigned int countInvalid = 0;
  for ( typename BcContainerType::const_iterator cit = m_vBc.begin();
        cit != m_vBc.end(); ++cit )
  {
    img = m_pmesh->dir_img( (*cit)->pt );
    if ( !img.isValid() )
    {
      ++countInvalid;
      continue;
    }
    // if a topology problem is observed, no point carrying on

    dval = ( img - (*cit)->pt - (*cit)->delta ).norm();

    dSumInitial += (*cit)->delta.norm();

    dSum += dval;
    ++count;
    if ( (*cit)->isActive )
    {
      dSumConditional += dval;
      countConditional++;
    }
  } // next cit
  std::cout << " countInvalid = " << countInvalid
  << " general-count = " << count << std::endl;

  if ( count )
  {
    std::cout << " Average of the error norm = "
    << dSum /(double)count << std::endl
    << " Initial error = " << dSumInitial / (double)count << std::endl;
  }
  else
    std::cout << " count = 0 !?!\n";

  if ( countConditional )
    std::cout << " Average of the error norm conditional = "
    << dSumConditional / (double)countConditional << std::endl;
  else
    std::cout << " countConditional = 0 !?!?\n";

  dRemainingRatio = dSum / dSumInitial;

  return 0;

}

//-----------------------------------------------------------------------
//
//

template<class Cstr, int n>
TDirectSolver<Cstr,n>::TDirectSolver()
    : TSolver<Cstr,n>()
{
  this->m_displayLevel = 1;
}

template<class Cstr, int n>
TDirectSolver<Cstr,n>::~TDirectSolver()
{}

#undef __FUNCT__
#define __FUNCT__ "TDirectSolver::solve"
template<class Cstr, int n>
int
TDirectSolver<Cstr,n>::solve()
{
  this->done_bc_natural();

  this->setup_matrix();

  // MFEM matrices are finalized during construction

  this->setup_load_sym();

  // MFEM vectors are ready to use after setup

  // Allocate solution vector
  this->m_delta = new mfem::Vector(this->m_load->Size());
  *(this->m_delta) = *(this->m_load); // Initialize with load vector

  // Create and configure the linear solver
  mfem::CGSolver cg;
  cg.SetRelTol(1.0e-9);
  cg.SetMaxIter(10000);
  cg.SetPrintLevel(this->m_displayLevel);
  cg.SetOperator(*(this->m_stiffness));

  // Solve the linear system
  cg.Mult(*(this->m_load), *(this->m_delta));

  if (cg.GetConverged())
  {
    if (this->m_displayLevel)
      std::cout << "DirectSolverConvergence = converged in "
                << cg.GetNumIterations() << " iterations" << std::endl;
  }
  else
  {
    std::cout << "WARNING: DirectSolver did not converge!" << std::endl;
  }

  // distribute obtained displacements to nodes
  this->comm_solution();

  // release MFEM resources
  delete this->m_delta;
  this->m_delta = nullptr;
  delete this->m_load;
  this->m_load = nullptr;
  delete this->m_stiffness;
  this->m_stiffness = nullptr;

  return 0;
}

#endif
