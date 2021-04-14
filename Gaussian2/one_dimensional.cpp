
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>

#include <meep.hpp>
#include <iostream>
#include <fstream>

using namespace meep;
using namespace std;
using std::complex;

typedef std::complex<double> cdouble;

double one(const vec &) { return 1.5; }
double two(const vec &) { return 0.1; }

cdouble line_integral_x(fields & f, component C, double dx, double xmin, double xmax, double y, double z)
{
  cdouble sum(0.0,0.0);
  cdouble deltax(dx,0.0);
  for (double x=xmin; x<=xmax; x=x+dx)
  {
    monitor_point p;
    f.get_point(&p, vec(x,y,z));
    cdouble dF = p.get_component(C);
    sum += dF*deltax;
  }
  return sum;
}

cdouble line_integral_y(fields & f, component C, double dy, double ymin, double ymax, double x, double z)
{
  cdouble sum(0.0,0.0);
  cdouble deltay(dy,0.0);
  for (double y=ymin; y<=ymax; y=y+dy)
  {
    monitor_point p;
    f.get_point(&p, vec(x,y,z));
    cdouble dF = p.get_component(C);
    sum += dF*deltay;
  }
  return sum;
}

cdouble line_integral_z(fields & f, component C, double dz, double zmin, double zmax, double x, double y)
{
  cdouble sum(0.0,0.0);
  cdouble deltaz(dz,0.0);
  for (double z=zmin; z<=zmax; z=z+dz)
  {
    monitor_point p;
    f.get_point(&p, vec(x,y,z));
    cdouble dF = p.get_component(C);
    //cout<<dF.real()<<" , "<<dF.imag()<<endl;
    sum += dF*deltaz;
    //cout<<dF.real()<<" , "<<dF.imag()<<endl;
  }
  return sum;
}



cdouble compute_Im(fields & f, double z)
{
    cdouble Iyf=line_integral_y(f,Ey,0.01,-1.0,1.0,1.0,z);
    cdouble Iyb=line_integral_y(f,Ey,0.01,-1.0,1.0,-1.0,z);
    cdouble Ixt=line_integral_x(f,Ex,0.01,-1.0,1.0,-1.0,z);
    cdouble Ixb=line_integral_x(f,Ex,0.01,-1.0,1.0,1.0,z);
    cdouble Im=Ixt+Iyf-Ixb-Iyb;
    return Im;
}


cdouble compute_Vm(fields & f, double z)
{
  //cdouble Vy=line_integral_y(f,Hy,0.001,ycen-dymin,y,xcen,zcen-dzmin);
  cdouble Vz=line_integral_y(f,Hy,0.01,-3.0,0.0,0.0,z);
  //cdouble Vm=Vy+Vz;
  return Vz;
}




int main(int argc, char **argv) {
  initialize mpi(argc, argv);
  verbosity = 4;
    std::ofstream FieldsIn;
    FieldsIn.open ("FieldEvolutionIn.txt");
    std::ofstream SourceFFT;
    SourceFFT.open ("SourceFFT.txt");
    std::ofstream Energy;
    Energy.open ("Energy.txt");
    std::ofstream Space;
    Space.open ("Space.txt");
    std::ofstream Voltages;
    Voltages.open ("Voltages.txt");

  double a = 24.0;
  double ttot = 100.0;

  grid_volume gv = vol3d(5.0,0.0,12.0, a);
  gv.center_origin();
  structure s(gv, one, pml(1.0), identity());
  s.add_susceptibility(two, H_stuff, gyrotropic_susceptibility(vec(0.0,0.0,1.0),1.0, 0.001,0.00001,GYROTROPIC_SATURATED));
  //bias omega gamma alpha type

  fields f(&s);
  f.use_real_fields();
  //continuous_src_time src(1.0);
  gaussian_src_time src(0.75,1.5);

for (double fp=0.0;fp<=2.0;fp=fp+0.0001)
    { 
      double yp = 0.0;
      {
        SourceFFT<<fp<<" "<<src.fourier_transform(fp).real()<<" "<<src.fourier_transform(fp).imag()<<endl;
      }
      
    }

  f.add_point_source(Hz,src,vec(0.0,0.0,-4.5));
    volume vin=volume(vec(-1,-1,-5),vec(1,1,-4));
    volume vout=volume(vec(-1,-1,4),vec(-1,-1,5));
    volume vyz=volume(vec(0,-6,-6),vec(0,6,6));
  int tind=0;
  while (f.time() < ttot) {
	f.step();tind++;
	//Energy<<f.magnetic_energy_in_box(vin)<<" , "<<f.magnetic_energy_in_box(vout)<<" , "<<f.magnetic_energy_in_box(gv.surroundings())<<endl;	
if((tind%1000)==0)
{
	for (double ind=-6.0;ind<=6.0;ind=ind+0.1)    
		{
    /*monitor_point pin;
        f.get_point(&pin, vec(0.0,0.0,ind));
        cdouble E1i = pin.get_component(Ex);
        cdouble E2i = pin.get_component(Ey);
        cdouble E3i = pin.get_component(Ez);
        cdouble D1i = pin.get_component(Dx);
        cdouble D2i = pin.get_component(Dy);
        cdouble D3i = pin.get_component(Dz);
        cdouble H1i = pin.get_component(Hx);
        cdouble H2i = pin.get_component(Hy);
        cdouble H3i = pin.get_component(Hz);
        cdouble B1i = pin.get_component(Bx);
        cdouble B2i = pin.get_component(By);
        cdouble B3i = pin.get_component(Bz);
        
        Space<<tind<<" , "<<ind<<" , "<<H1i.real() <<" , "<<H1i.imag()<<" , "<<H2i.real()<<" , "<<H2i.imag()<<" , "<<H3i.real()<<" , "<<H3i.imag()<<" , "<<B1i.real()<<" , "<<B1i.imag()<<" , "<<B2i.real()<<" , "<<B2i.imag()<<" , "<<B3i.real()<<" , "<<B3i.imag()<<" , "<<E1i.real()<<" , "<<E1i.imag()<<" , "<<E2i.real()<<" , "<<E2i.imag()<<" , "<<E3i.real()<<" , "<<E3i.imag()<<" , "<<D1i.real()<<" , "<<D1i.imag()<<" , "<<D2i.real()<<" , "<<D2i.imag()<<" , "<<D3i.real()<<" , "<<D3i.imag()<<" , "<<endl;		

        cdouble Vmi=compute_Vm(f,ind);
        cdouble Imi=compute_Im(f,ind);
        Voltages<<tind<<" , "<<ind<<" , "<<Imi.real()<<" , "<<Imi.imag()<<" , "<<Vmi.real()<<" , "<<Vmi.imag()<<endl;
//" , "<<Ie.real()<<" , "<<Ie.imag()<<" , "<<Ve.real()<<" , "<<Ve.imag()<<endl;
*/
}
}



	for (double ind=-6.0;ind<=6.0;ind=ind+1)    
		{
    monitor_point pin;
        f.get_point(&pin, vec(0.0,0.0,ind));
        cdouble E1i = pin.get_component(Ex);
        cdouble E2i = pin.get_component(Ey);
        cdouble E3i = pin.get_component(Ez);
        cdouble D1i = pin.get_component(Dx);
        cdouble D2i = pin.get_component(Dy);
        cdouble D3i = pin.get_component(Dz);
        cdouble H1i = pin.get_component(Hx);
        cdouble H2i = pin.get_component(Hy);
        cdouble H3i = pin.get_component(Hz);
        cdouble B1i = pin.get_component(Bx);
        cdouble B2i = pin.get_component(By);
        cdouble B3i = pin.get_component(Bz);
        
        FieldsIn<<tind<<" , "<<ind<<" , "<<H1i.real() <<" , "<<H1i.imag()<<" , "<<H2i.real()<<" , "<<H2i.imag()<<" , "<<H3i.real()<<" , "<<H3i.imag()<<" , "<<B1i.real()<<" , "<<B1i.imag()<<" , "<<B2i.real()<<" , "<<B2i.imag()<<" , "<<B3i.real()<<" , "<<B3i.imag()<<" , "<<E1i.real()<<" , "<<E1i.imag()<<" , "<<E2i.real()<<" , "<<E2i.imag()<<" , "<<E3i.real()<<" , "<<E3i.imag()<<" , "<<D1i.real()<<" , "<<D1i.imag()<<" , "<<D2i.real()<<" , "<<D2i.imag()<<" , "<<D3i.real()<<" , "<<D3i.imag()<<" , "<<endl;		
}
  //complex<double> E1i = pin.get_component(Ez);  
  //cout<<E1i.real()<<endl;
//f.output_hdf5(Hy,vxz);
//f.output_hdf5(Hx,vxz);

	}
  return 0;
}
