#include "ReadFromXML.h"
#include "Geometry.h"
#include "log.h"
#include "CPUSolver.h"
#include "TrackGenerator3D.h"
#include "Quadrature.h"
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

using namespace antmoc;

int main(){
    Geometry *g = new Geometry();
    ReadFromXML *inp =new ReadFromXML(g);
    inp->readMaterials("cases/Beavrs/mac.h5");
    g=inp->readGeometry("cases/Beavrs/Beavrs_3D.xml");
    Universe *root_universe = g->getRootUniverse();
    //std::cout<<root_universe->toString()<<std::endl;
	//构建完了几何就对几何进行扇区的划分
    log_printf(NORMAL,"devide the cells");
    g->initializeFlatSourceRegions();
	log_printf(NORMAL,"set material to the Material Cells");
	inp->setMaterialsToCells();

    //设置参数
     #ifdef OPENMP
	 int num_threads = omp_get_num_procs();
	  #else
	  int num_threads = 1;
      #endif
	  double azim_spacing = 0.5;
	  int num_azim = 4;
	  double polar_spacing = 0.4;
	  int num_polar = 4;
	  double tolerance = 1e-5;
	  int max_iters = 30;
	  int axial_refines = 1;
	  log_printf(NORMAL, "Initializing the track generator...");
	 // Quadrature* quad = new GLPolarQuad();
	 Quadrature* quad = new EqualAnglePolarQuad();
	  quad->setNumPolarAngles(num_polar);
	 // quad->setNumAzimAngles(num_azim);
	  TrackGenerator3D track_generator(g, num_azim, num_polar, azim_spacing,
	                                      polar_spacing);
	  track_generator.setNumThreads(num_threads);
	  track_generator.setQuadrature(quad);
	  track_generator.setSegmentFormation(OTF_STACKS);
	  std::vector<double> seg_zones{-500,500};
	  track_generator.setSegmentationZones(seg_zones);
	  track_generator.generateTracks();


	   log_printf(NORMAL, "run simulation...");
	   CPUSolver solver(&track_generator);
	   solver.setNumThreads(num_threads);
	   solver.setConvergenceThreshold(tolerance);
	   solver.computeEigenvalue(max_iters);
	    solver.printTimerReport();
		log_printf(TITLE, "Finished");

}

