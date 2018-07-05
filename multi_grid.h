/* ----------------------------------------------------------------------
Poisson solver subclass
------------------------------------------------------------------------- */
#ifndef LMP_MULTI_GRID
#define LMP_MULTI_GRID

#include <complex>
#include <random>

namespace LAMMPS_NS {

class MULTI_GRID{
    public:
    
    MULTI_GRID(int *, double *, double *, double, bool *);
    ~MULTI_GRID();
    
    double *** array, *** array_g;
    double * gs, * lo, *hi, *gis;
    int * nt;
    double lambda, sigma;
    
    void gradient_calc(double *, double *);
    
    void ens_average(MPI_Comm *);
    //void update_boundary();
    
    void zero();
    void subtract_av();
    void scale(double);
    
    double lap_calc(int, int, int);
    
    
    void write_slab(double ***, FILE *, bigint, int);
    void write_rad(double ***, FILE *, bigint, int,double,double);
    
    bool ind_calc(int *, double *);
    void ind_calc_floor(int *, double *);
    void pos_calc(int *, double *);
    
    void add_sc_up(MULTI_GRID *);
    void add_sc_down(MULTI_GRID *);
    void spread_charge(double);
    
    void update_boundary(int, int, int);
    void heat_iterate(double *** , double, double, bool);
        
    private:
    
    void gauss_conv(int dim, int na, int nb, double st);
    int return_opposite(int, int);
    void delete_array(int *, double ***);
    bool slab_flag, ens, offset_flag;
    double *** create_array(int *);
    int mod(int, int);
    bool div_check(int, int);
    double compute_my_gfunc(double);
    double calculate_pf(int, double, double, double);
    
};

}

#endif