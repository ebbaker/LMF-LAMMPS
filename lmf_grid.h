/* ----------------------------------------------------------------------
Poisson solver subclass
------------------------------------------------------------------------- */
#ifndef LMP_LMF_GRID
#define LMP_LMF_GRID

#include <complex>

namespace LAMMPS_NS {

class LMF_GRID{
    public:
    
    LMF_GRID(int *, double *, double *, bool);
    ~LMF_GRID();
    
    double *** array, *** array_g;
    
    int *nt;
    double *gs, *gis, *lo;
    
    void gradient_calc(int *, double *);
    void exp_smooth(int *, int *, double, double, double, double *, double *);
    void ens_average(MPI_Comm *);
	void sc_update();
    
    void zero();
    void subtract_av(double);
    
    double lap_calc(int , int, int);
    
    private:
    
    bool ens;
    double *** create_array(int *);
};

}

#endif