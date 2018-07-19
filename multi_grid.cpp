/* ----------------------------------------------------------------------
Poisson solver subclass
------------------------------------------------------------------------- */
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "variable.h"
#include "multi_grid.h"
#include <complex>
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "atom.h"
#include "input.h"
#include "variable.h"
#include "domain.h"
#include "lattice.h"
#include "update.h"
#include "modify.h"
#include "respa.h"
#include "error.h"
#include "force.h"
#include "universe.h"
#include <complex>
#include "run.h"
#include <random>

using namespace LAMMPS_NS;
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

MULTI_GRID::MULTI_GRID(int * nt_t, double * gs_t, double * lo_t, double sigma_t, bool * flags){
    ens=flags[0]; slab_flag=flags[1]; 
    sigma=sigma_t;
    nt=nt_t; gs=gs_t; lo=lo_t;  
    array=create_array(nt);
    if(ens) array_g=create_array(nt);
    else array_g=NULL;
    gis=new double[3]; hi=new double[3];
    for(int dim=0; dim<3; dim++){
        gis[dim]=1.0/gs[dim]/gs[dim];
        hi[dim]=lo[dim]+double(nt[dim]-1)*gs[dim];
    }
}


MULTI_GRID::~MULTI_GRID(){
    delete_array(nt,array);
    if(ens)delete_array(nt,array_g);
    delete [] gis; delete [] gs; 
    delete [] lo; delete [] hi;
    delete [] nt;
}

void MULTI_GRID::gradient_calc(double * x, double * gr){
    double base, ltemp=0;
    int ind[3]; 
    ind_calc_floor(ind,x);
    double base_potential = array[ind[0]][ind[1]][ind[2]];
    gr[0]=(array[ind[0]+1][ind[1]][ind[2]]-base_potential)/gs[0];
	gr[1]=(array[ind[0]][ind[1]+1][ind[2]]-base_potential)/gs[1];
	gr[2]=(array[ind[0]][ind[1]][ind[2]+1]-base_potential)/gs[2];
}


void MULTI_GRID::ens_average(MPI_Comm *uorig){
    if(ens){
        int num_procs;
        double temp_array[nt[2]];
        MPI_Comm_size (*uorig, &num_procs);
        for(int n_0=0; n_0<nt[0]; n_0++){
            for(int n_1=0; n_1<nt[1]; n_1++){
                MPI_Allreduce(array[n_0][n_1],temp_array,nt[2],MPI_DOUBLE,MPI_SUM,*uorig);
                for(int n_2=0; n_2<nt[2]; n_2++){
                    array[n_0][n_1][n_2]=temp_array[n_2]/double(num_procs);
                }
            }
        }
    }
}


void MULTI_GRID::zero(){
    for(int n_0=0; n_0<nt[0]; n_0++){
        for(int n_1=0; n_1<nt[1]; n_1++){
            for(int n_2=0; n_2<nt[2]; n_2++){
                array[n_0][n_1][n_2]=0.0;
                if (ens) array_g[n_0][n_1][n_2]=0.0;
            }	
        }		
    }	
}


bool MULTI_GRID::ind_calc(int * ind, double * x){
    bool return_val=1;
    for(int dim=0; dim<3; dim++){
        ind[dim] = int(round((x[dim]-lo[dim])/gs[dim]))+1;
        if(!div_check(ind[dim],nt[dim]-1)){
            return_val=0;
        }
        else if(ind[dim]<0 || ind[dim]>=nt[dim]-1){
            if(dim==2 && slab_flag) return_val=0;   
            ind[dim]=mod(ind[dim],nt[dim]-2);
        }
    }
    return return_val;
}

void MULTI_GRID::ind_calc_floor(int * ind, double * x){
    for(int dim=0; dim<3; dim++){
        ind[dim] = int(floor((x[dim]-lo[dim])/gs[dim]))+1;
        if(ind[dim]<0 || ind[dim]>=nt[dim]-1){
            ind[dim]=mod(ind[dim],nt[dim]-2);
        }
    }
}

void MULTI_GRID::add_sc_up(MULTI_GRID * sc){
    int i[3]; double x[3]; double xf[3]; double base;
    double pf[3]; double gs_r=1.0;
    for(int dim=0; dim<3; dim++) gs_r*=float(sc->nt[dim]-2)/float(nt[dim]-2)/sc->gs[dim];
    for(int n0=1; n0<nt[0]-1; n0++){
		for(int n1=1; n1<nt[1]-1; n1++){
			for(int n2=1; n2<nt[2]-1; n2++){
			    i[0]=n0; i[1]=n1; i[2]=n2;
			    pos_calc(i,x);
			    sc->ind_calc_floor(i, x);
			    sc->pos_calc(i,xf);
			    base=array[n0][n1][n2];
			    for(int w0=0;w0<2;w0++){
			        pf[0]=calculate_pf(w0, x[0], xf[0], sc->gs[0]);
			        for(int w1=0;w1<2;w1++){
			            pf[1]=calculate_pf(w1, x[1], xf[1], sc->gs[1]);
			            for(int w2=0;w2<2;w2++){
			                pf[2]=calculate_pf(w2, x[2], xf[2], sc->gs[2]);
			                sc->array[i[0]+w0][i[1]+w1][i[2]+w2]+=base*pf[0]*pf[1]*pf[2]*gs_r;
			                sc->update_boundary(i[0]+w0,i[1]+w1,i[2]+w2);
			                //printf("%d %d %d %d %d %d %d %d %d %f \n", n0, n1, n2, i[0], i[1], i[2], w0, w1, w2, pf[0]*pf[1]*pf[2]*gs_r);
			            }
			        }
			    }
			}
	    }
	}   
}

void MULTI_GRID::add_sc_down(MULTI_GRID * sc){
    int i[3]; double x[3]; double xf[3]; double base;
    int w[3]; double pf[3]; double gs_i=1.0;
    for(int dim=0; dim<3; dim++) gs_i/=sc->gs[dim];
    for(int n0=1; n0<nt[0]-1; n0++)
		for(int n1=1; n1<nt[1]-1; n1++)
			for(int n2=1; n2<nt[2]-1; n2++){
			    i[0]=n0; i[1]=n1; i[2]=n2;
			    pos_calc(i,x);
			    sc->ind_calc_floor(i, x);
			    sc->pos_calc(i,xf);
			    for(int w0=0;w0<2;w0++){
			        pf[0]=calculate_pf(w0, x[0], xf[0], sc->gs[0]);
			        for(int w1=0;w1<2;w1++){
			            pf[1]=calculate_pf(w1, x[1], xf[1], sc->gs[1]);
			            for(int w2=0;w2<2;w2++){
			                pf[2]=calculate_pf(w2, x[2], xf[2], sc->gs[2]);
			                array[n0][n1][n2]-=sc->array[i[0]+w0][i[1]+w1][i[2]+w2]*pf[0]*pf[1]*pf[2]*gs_i;
			                update_boundary(n0,n1,n2);
			            }
			        }
			    }
			}
}

double MULTI_GRID::calculate_pf(int w, double x, double xf, double gs){
    if(abs(x-xf)>gs) printf("Error with calculation %d %f %f %f\n", w, x, xf, gs);
    else if(w==0) return (gs-x+xf);
    else if(w==1) return (x-xf);
}

void MULTI_GRID::spread_charge(double st){
    
    for(int n_a=1; n_a<nt[1]-1; n_a++){
        for(int n_b=1; n_b<nt[2]-1; n_b++){
            gauss_conv(0, n_a, n_b, st);
        }
    }
    
    for(int n_a=1; n_a<nt[0]-1; n_a++){
        for(int n_b=1; n_b<nt[2]-1; n_b++){
            gauss_conv(1, n_a, n_b, st);
        }
    }
    
	for(int n_a=1; n_a<nt[0]-1; n_a++){
        for(int n_b=1; n_b<nt[1]-1; n_b++){
            gauss_conv(2, n_a, n_b, st);
        }
    }
	
}

void MULTI_GRID::heat_iterate(double *** sct, double hpft, double lpf, bool error_calc){
    double total_error=0.0;
    for(int n_0=1; n_0<nt[0]-1; n_0++){
        for(int n_1=1; n_1<nt[1]-1; n_1++){
            for(int n_2=1; n_2<nt[2]-1; n_2++){
                array[n_0][n_1][n_2]+=hpft*(lap_calc(n_0,n_1,n_2)+lpf*sct[n_0][n_1][n_2]);
                update_boundary(n_0,n_1,n_2);
                if(error_calc){
                    double error=(lap_calc(n_0,n_1,n_2)+lpf*sct[n_0][n_1][n_2])*(lap_calc(n_0,n_1,n_2)+lpf*sct[n_0][n_1][n_2]);
                    total_error+=error;
                    if(error>10.0) printf("error %d %d %d %f %f %f %e \n",n_0,n_1,n_2, array[n_0][n_1][n_2], lap_calc(n_0,n_1,n_2), lpf*sct[n_0][n_1][n_2], error);
                }
            }
        }
    }
    if(error_calc) printf("error %e \n", total_error);
	subtract_av();
} 


void MULTI_GRID::gauss_conv(int dim, int na, int nb, double st){
    double gst = gs[dim]; int nm = nt[dim]; 
    int max = int(3.0*st/gst)+1; int ind_t;
    double at[nm], gauss_t[2*max+1], conv_t[nm];
    
    double norm=0.0;
    for(int n=-max; n<max; n++){
        gauss_t[max+n]=exp(-double(n*n)*gst*gst/st/2.0);
        norm+=gauss_t[max+n];
    }
    for(int n=-max; n<max+1; n++) gauss_t[max+n]=gauss_t[max+n]/norm;

    if(dim==0) for(int i=0; i<nm; i++) at[i]=array[i][na][nb];
    if(dim==1) for(int i=0; i<nm; i++) at[i]=array[na][i][nb];
    if(dim==2) for(int i=0; i<nm; i++) at[i]=array[na][nb][i];
        
    for(int i=1; i<nm-1; i++){    
        conv_t[i]=0.0;
        for(int j=i-max; j<i+max+1; j++){
            ind_t = mod(j,nm-2);
            conv_t[i]+=at[ind_t]*gauss_t[j-i+max];
        }
    }
    conv_t[nm-1]=conv_t[1]; conv_t[0]=conv_t[nm-2];
    
    if(dim==0) for(int i=0; i<nm; i++){
                    array[i][na][nb]=conv_t[i];
                    update_boundary(i,na,nb);
                }
    if(dim==1) for(int i=0; i<nm; i++){
                    array[na][i][nb]=conv_t[i];
                    update_boundary(na,i,nb);
                }
    if(dim==2) for(int i=0; i<nm; i++){
                    array[na][nb][i]=conv_t[i];
                    update_boundary(na,nb,i);
                }
} 

void MULTI_GRID::pos_calc(int * ind, double * x){
    for(int dim=0; dim<3; dim++){
        x[dim]=double(ind[dim]-1)*gs[dim]+lo[dim];
    }
} 

int MULTI_GRID::mod(int a, int b)
{
    if(a<0) a=a-(a/b-1)*b;
    return(a%b);
}

bool MULTI_GRID::div_check(int a, int b)
{
   if (a>2*b) return 0;
   else if(-a>b) return 0;
   else return 1;
}


double MULTI_GRID::lap_calc(int n_0, int n_1, int n_2){
    double temp=0.0;
    double baset = 2.0*array[n_0][n_1][n_2];
	temp += (array[n_0+1][n_1][n_2]+array[n_0-1][n_1][n_2]-baset)*gis[0];
	temp += (array[n_0][n_1+1][n_2]+array[n_0][n_1-1][n_2]-baset)*gis[1];
	temp += (array[n_0][n_1][n_2+1]+array[n_0][n_1][n_2-1]-baset)*gis[2];
	return temp;
}


double *** MULTI_GRID::create_array(int * size){
	double *** ar = new double ** [size[0]];
	for(int n0=0; n0<size[0]; n0++){
		ar[n0]=new double * [size[1]];
		for(int n1=0; n1<size[1]; n1++){
			ar[n0][n1]= new double [size[2]];
			for(int n2=0; n2<size[2]; n2++)
				ar[n0][n1][n2]=0.0;
		}
	}
	return ar;
}

void MULTI_GRID::delete_array(int * size, double *** array){
	for(int n0=0; n0<size[0]; n0++){
		for(int n1=0; n1<size[1]; n1++)
			delete [] array[n0][n1];
		delete [] array[n0];
	}
	delete [] array;
}

void MULTI_GRID::subtract_av(){
    double total_field=0.0;
    for(int n_0=0; n_0<nt[0]; n_0++){
		for(int n_1=0; n_1<nt[1]; n_1++){
			for(int n_2=0; n_2<nt[2]; n_2++){
				total_field+=array[n_0][n_1][n_2];
			}	
		}		
	}	
	total_field=total_field/double(nt[0]*nt[1]*nt[2]);
	for(int n_0=0; n_0<nt[0]; n_0++){
		for(int n_1=0; n_1<nt[1]; n_1++){
			for(int n_2=0; n_2<nt[2]; n_2++){
				array[n_0][n_1][n_2]-=total_field;
			}	
		}		
	}		
}

void MULTI_GRID::scale(double scale_val){
    for(int n_0=0; n_0<nt[0]; n_0++){
		for(int n_1=0; n_1<nt[1]; n_1++){
			for(int n_2=0; n_2<nt[2]; n_2++){
				array[n_0][n_1][n_2]*=scale_val;
			}	
		}		
	}
}

void MULTI_GRID::write_slab(double *** array_2, FILE * ft, bigint ts, int gn){
    double temp_val, temp_val_2;
    fprintf(ft,"Timestep %d %d %d %f %f %f \n",ts,gn,nt[2],gs[2],lo[2],hi[2]);
    for(int n_2=0; n_2<nt[2]; n_2++){
        temp_val=0.0; temp_val_2=0.0;
        for(int n_0=1; n_0<nt[0]-1; n_0++){
	    	for(int n_1=1; n_1<nt[1]-1; n_1++){
		    	temp_val+=array[n_0][n_1][n_2];
		    	temp_val_2+=array_2[n_0][n_1][n_2];
			}	
		}	
		temp_val=temp_val/double((nt[0]-2)*(nt[1]-2));
		temp_val_2=temp_val_2/double((nt[0]-2)*(nt[1]-2));
		fprintf(ft,"%d %f %f \n",n_2,temp_val,temp_val_2);
	}
	fflush(ft);
}


void MULTI_GRID::write_rad(double *** array_2, FILE * ft, bigint ts, int gn, double max_rad, double spac){
    double temp_val, temp_val_2;
    int max_ind = int(max_rad/spac);
    printf("gn %d max_ind %d max_rad %f spac %f \n",gn, max_ind, max_rad, spac);
    fprintf(ft,"Timestep %d %d %d %f \n",ts,gn,max_ind,spac);
    //fprintf(ft,"Timestep %d %d %d %f \n",ts,gn,nt[0],gs[0]);
    int indt[3]; double xt[3];
    int count_array[max_ind];
    double sc_vals[max_ind];
    double pot_vals[max_ind];
    double lmf_pf = 180.951271/4.0/3.1415;
    for(int n_r=0; n_r<max_ind; n_r++){
        count_array[n_r]=0;
        sc_vals[n_r]=0.0;
        pot_vals[n_r]=0.0;
    }
    for(int n_0=1; n_0<max_rad/gs[0]; n_0++)
        for(int n_1=1; n_1<max_rad/gs[1]; n_1++)
            for(int n_2=1; n_2<max_rad/gs[2]; n_2++){
                indt[0]=n_0; indt[1]=n_1; indt[2]=n_2;
                pos_calc(indt,xt);
                double dist = sqrt(pow(xt[0],2.0)+pow(xt[1],2.0)+pow(xt[2],2.0));
                int dind = int(dist/spac);
                if(dind<max_ind){
                    sc_vals[dind]=sc_vals[dind]+array[n_0][n_1][n_2];
                    pot_vals[dind]=pot_vals[dind]+array_2[n_0][n_1][n_2];
                    count_array[dind]=count_array[dind]+1;
                }
            }
    for(int n_r=0; n_r<max_ind; n_r++){
        if(count_array[n_r]!=0){
            sc_vals[n_r]=sc_vals[n_r]/float(count_array[n_r]);
            pot_vals[n_r]=pot_vals[n_r]/float(count_array[n_r]);
        }
		fprintf(ft,"%d %f %f %f \n",n_r,float(n_r)*spac, sc_vals[n_r],pot_vals[n_r]);
	}
	//for(int n_0=0; n_0<nt[0]; n_0++)
	//    fprintf(ft,"%d %f %f %f \n",n_0,float(n_0)*gs[0], array[n_0][1][1], array_2[n_0][1][1]);
	fflush(ft);
}


void MULTI_GRID::update_boundary(int n_0, int n_1, int n_2){
    int opp_0, opp_1, opp_2;
    double base = array[n_0][n_1][n_2];
    opp_0=return_opposite(n_0,nt[0]);
    opp_1=return_opposite(n_1,nt[1]);
    opp_2=return_opposite(n_2,nt[2]);
    if(opp_0!=n_0) array[opp_0][n_1][n_2]=base;
    if(opp_1!=n_1){
        array[n_0][opp_1][n_2]=base;
        if(opp_0!=n_0) array[opp_0][opp_1][n_2]=base;
        }
    if(opp_2!=n_2){
        if(slab_flag){
            if(n_2==1) opp_2=0;
            if(n_2==nt[2]-2) opp_2=nt[2]-1;
        }
        array[n_0][n_1][opp_2]=base; 
        if(opp_0!=n_0) array[opp_0][n_1][opp_2]=base; 
        if(opp_1!=n_1) array[n_0][opp_1][opp_2]=base;
        if(opp_0!=n_0 && opp_1!=n_1) array[opp_0][opp_1][opp_2]=base;
    }
}

int MULTI_GRID::return_opposite(int ind,int high){
    if(ind>1 && ind<high-2) return ind;
    else if(ind==1) return (high-1);
    else if(ind==(high-2)) return 0;
    else if(ind==0) return (high-2);
    else if(ind==(high-1)) return 1;
	else return -1;
}



