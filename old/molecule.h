
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"

using namespace std;

class Molecule {
  public:
    int natom;
    int charge;
    int *zvals;
    double **geom;
    double atomic_mass[];

    void print_geom();
    void rotate(double phi);
    void translate(double x, double y, double z);
    void set_geom(const char *filename);
    
    double unit(int cart, int atom1, int atom2);
    void unitvec(int atom1, int atom2, Eigen::Vector3d v);
    double bond(int atom1, int atom2);
    double angle(int atom1, int atom2, int atom3);
    double outofplane(int atom1, int atom2, int atom3, int atom4);
    double torsional(int atom1, int atom2, int atom3, int atom4);
    double masscenter(int cart);
    void inertia(Eigen::Matrix3d *evals, Eigen::Matrix3d *evecs);

    Molecule(const char *filename, int q);
    ~Molecule();
};
