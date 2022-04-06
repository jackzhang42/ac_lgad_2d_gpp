#include <iostream>
#include <fstream>
#include <cmath>

#include <string>
#include <vector>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>

#include <TTree.h>
#include <TFile.h>

#include "Garfield/MediumSilicon.hh"
#include "Garfield/ComponentTcad2d.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewSignal.hh"
#include "Garfield/FundamentalConstants.hh"
#include "Garfield/Random.hh"

using namespace Garfield;

// define transfer function
double transfer(double t) {
	constexpr double tauR = 5.6;
	constexpr double tauI = 1.8;
	constexpr double tauA = 47.;
	constexpr double dAI = (tauA - tauI);
	constexpr double dAR = (tauA - tauR);
	constexpr double dIR = (tauI - tauR);
	constexpr double c1 = tauA / (dAI * dAI *dAR);
	constexpr double c2 = 1. / (dAI * tauI * dIR);
	constexpr double c3 = tauR / (dAR * dIR * dIR);
	constexpr double c4 = (tauI * tauI * tauI - tauA * tauI * tauR) / (dAI * dAI * tauI * dIR * dIR);
	const double f1 = -exp(-t / tauA) * c1;
	const double f2 = exp(-t / tauI) * t * c2;
	const double f3 = exp(-t / tauR) * c3;
	const double f4 = exp(-t / tauI) * c4;
	return tauA * tauR * (f1 + f2 + f3 + f4);

}

int main(int argc, char * argv[]) {
	
	TApplication app("app", &argc, argv);

// define constants
	const double depth = 50.e-4;
	const double left = 175.e-4;
	const double right = 1125.e-4;
	const double width = 650.e-4;
	const double yhigh = 1300.e-4;
  const int npad = 5;
	const double gap = 75.e-4;
	const double side = 100.e-4;
	const double tmetal = 1.e-4;
	const double lnp = npad * side + (npad + 1) * gap + 2 * tmetal;
	const double npleft = (yhigh - lnp) / 2;
	const double npright = (yhigh + lnp) / 2;
	const double acleft = npleft + tmetal;
	const double acright = npright - tmetal;
	const double pad1left = acleft + gap;
	const double pad1right = pad1left + side;
	const double pad2left = pad1right + gap;
	const double pad2right = pad2left + side;
	const double pad3left = pad2right + gap;
	const double pad3right = pad3left + side;
	const double pad4left = pad3right + gap;
	const double pad4right = pad4left + side;
	const double pad5left = pad4right + gap;
	const double pad5right = pad5left + side;

// define flags
	constexpr bool plotField = false;
  constexpr bool plotWeightingField = false;
	constexpr bool plotSignal = true;
	constexpr bool plotDrift = false;
  constexpr bool writeSignal = true;
	constexpr bool writeTime = false;
	constexpr bool useTransfer = false;

	MediumSilicon si;
	ComponentTcad2d cmp;
//  cmp.EnableDebugging(); // open debug mode
//  cmp.Initialise("normal_bias_3253_0000_lgad_des.grd", "normal_bias_3253_0000_lgad_des.dat");
  cmp.Initialise("normal_bias_10_0000_lgad_des.grd", "normal_bias_10_0000_lgad_des.dat");


// weighting field from tcad file
	ComponentTcad2d dwp_1;
	ComponentTcad2d dwp_2;
	ComponentTcad2d dwp_3;
	ComponentTcad2d dwp_4;
	ComponentTcad2d dwp_5;
	double dv = 1;
	dwp_1.Initialise("normal_bias_10_0000_lgad_des.grd", "normal_bias_10_0000_lgad_des.dat");
	dwp_2.Initialise("normal_bias_10_0000_lgad_des.grd", "normal_bias_10_0000_lgad_des.dat");
	dwp_3.Initialise("normal_bias_10_0000_lgad_des.grd", "normal_bias_10_0000_lgad_des.dat");
	dwp_4.Initialise("normal_bias_10_0000_lgad_des.grd", "normal_bias_10_0000_lgad_des.dat");
	dwp_5.Initialise("normal_bias_10_0000_lgad_des.grd", "normal_bias_10_0000_lgad_des.dat");
	dwp_1.SetWeightingField("normal_bias_10_0000_lgad_des.dat", "dynamic_n10_0000_lgad_des.dat", dv, "pad1-delayed");
	dwp_2.SetWeightingField("normal_bias_10_0000_lgad_des.dat", "dynamic_n11_0000_lgad_des.dat", dv, "pad2-delayed");
	dwp_3.SetWeightingField("normal_bias_10_0000_lgad_des.dat", "dynamic_n12_0000_lgad_des.dat", dv, "pad3-delayed");
	dwp_4.SetWeightingField("normal_bias_10_0000_lgad_des.dat", "dynamic_n13_0000_lgad_des.dat", dv, "pad4-delayed");
	dwp_5.SetWeightingField("normal_bias_10_0000_lgad_des.dat", "dynamic_n14_0000_lgad_des.dat", dv, "pad5-delayed");
	
	double start = 4.e-3;
	double end = 5.;
	double q = pow(end / start, 1./50);
	std::vector<double> dy_times;
	std::string prefix_1 = "dynamic_n10_";
	std::string prefix_2 = "dynamic_n11_";
	std::string prefix_3 = "dynamic_n12_";
	std::string prefix_4 = "dynamic_n13_";
	std::string prefix_5 = "dynamic_n14_";
	std::string profix = "_lgad_des.dat";

	for (int i = 0; i < 51; i++) {
		double t = start*pow(q, i);
		dy_times.push_back(t);
		std::string filename;
		if (i < 10) filename = prefix_1 + "000" + std::to_string(i) + profix;
		if (i > 9 && i < 100) filename = prefix_1 + "00" + std::to_string(i) + profix;
		dwp_1.SetWeightingField("normal_bias_10_0000_lgad_des.dat", filename, dv, t, "pad1-delayed");
		std::cout << "pad1: " << "t = " << t << ", happy" << i << "\n";
	}
	for (int i = 0; i < 51; i++) {
		double t = start*pow(q, i);
		std::string filename;
		if (i < 10) filename = prefix_2 + "000" + std::to_string(i) + profix;
		if (i > 9 && i < 100) filename = prefix_2 + "00" + std::to_string(i) + profix;
		dwp_2.SetWeightingField("normal_bias_10_0000_lgad_des.dat", filename, dv, t, "pad2-delayed");
		std::cout << "pad2: " << "t = " << t << ", happy" << i << "\n";
	}
	for (int i = 0; i < 51; i++) {
		double t = start*pow(q, i);
		std::string filename;
		if (i < 10) filename = prefix_3 + "000" + std::to_string(i) + profix;
		if (i > 9 && i < 100) filename = prefix_3 + "00" + std::to_string(i) + profix;
		dwp_3.SetWeightingField("normal_bias_10_0000_lgad_des.dat", filename, dv, t, "pad3-delayed");
		std::cout << "pad3: " << "t = " << t << ", happy" << i << "\n";
	}
	for (int i = 0; i < 51; i++) {
		double t = start*pow(q, i);
		std::string filename;
		if (i < 10) filename = prefix_4 + "000" + std::to_string(i) + profix;
		if (i > 9 && i < 100) filename = prefix_4 + "00" + std::to_string(i) + profix;
		dwp_4.SetWeightingField("normal_bias_10_0000_lgad_des.dat", filename, dv, t, "pad4-delayed");
		std::cout << "pad4: " << "t = " << t << ", happy" << i << "\n";
	}
	for (int i = 0; i < 51; i++) {
		double t = start*pow(q, i);
		std::string filename;
		if (i < 10) filename = prefix_5 + "000" + std::to_string(i) + profix;
		if (i > 9 && i < 100) filename = prefix_5 + "00" + std::to_string(i) + profix;
		dwp_5.SetWeightingField("normal_bias_10_0000_lgad_des.dat", filename, dv, t, "pad5-delayed");
		std::cout << "pad5: " << "t = " << t << ", happy" << i << "\n";
	}


	cmp.SetRangeZ(-width, width);

// assign silicon to corresponding region
	int nRegions = cmp.GetNumberOfRegions();
	for (int i = 0; i < nRegions; ++i) {
		std::string region;
		bool active;
		cmp.GetRegion(i, region, active);
		if (region == "\"Silicon_1\"") cmp.SetMedium(i, &si);
	}

// view potential field	
	ViewField vField;
	if (plotField) {
		vField.SetComponent(&cmp);
		vField.PlotContour("v");
	}

// view weightfield	
  if (plotWeightingField) {
    ViewField* wfieldView = new ViewField();
		wfieldView->SetComponent(&dwp_1);
    wfieldView->SetArea(-2e-4, 0, depth, 1300.e-4);
		wfieldView->PlotContourWeightingField("pad1-delayed", "v");
	}

	Sensor sensor;

//	sensor.EnableDebugging();

	sensor.AddComponent(&cmp);
	
	sensor.AddElectrode(&dwp_1, "pad1-delayed");
	sensor.AddElectrode(&dwp_2, "pad2-delayed");
	sensor.AddElectrode(&dwp_3, "pad3-delayed");
	sensor.AddElectrode(&dwp_4, "pad4-delayed");
	sensor.AddElectrode(&dwp_5, "pad5-delayed");

	const std::string label1 = "strip1";
	const std::string label2 = "strip2";
	const std::string label3 = "strip3";
	const std::string label4 = "strip4";
	const std::string label5 = "strip5";
	const std::string label6 = "cathode";
	const std::string label7 = "cathode1";

	double endtime = 5;
//	double endtime = 4.e-3;
	if (useTransfer) {
		endtime = 20;
		sensor.SetTransferFunction(transfer);
	}
//	const int nSignalBins = 100000;
	const int nSignalBins = 100;
	double tStep = endtime / nSignalBins;
	sensor.SetTimeWindow(0., tStep, nSignalBins);
	
	sensor.EnableDelayedSignal();
	
	std::vector<double> times;
	for (unsigned int i = 0; i < nSignalBins; i++) {
		times.push_back(tStep / 2 + i * tStep);
	}
  sensor.SetDelayedSignalTimes(times);
	const double thr1 = -1000 * ElementaryCharge;
	std::cout << "Threshold: " << thr1 << " fC\n";

	AvalancheMC drift;
//	drift.SetDistanceSteps(0.1e-4);
//	drift.SetDistanceSteps(1e-4);
//	drift.SetTimeSteps(0.05);
	drift.SetTimeSteps(0.1);
	drift.SetSensor(&sensor);
	drift.EnableSignalCalculation();
	drift.UseWeightingPotential();
//	drift.EnableInducedChargeCalculation();

	TrackHeed track;
	track.SetSensor(&sensor);
	track.SetParticle("pi");
	track.SetMomentum(420.e6); // MIP

// view signal
	ViewSignal vSignal;
	vSignal.SetSensor(&sensor);
	
// view drift lines	
	ViewDrift vDrift;
	if (plotDrift) {
		vDrift.SetArea(0., 0., depth, 1300.e-4);
		track.EnablePlotting(&vDrift);
	}

// open files	
	std::ofstream output_charge, output_time, output_gain;
	std::ofstream log;
	std::ofstream output_signal, output_y_real;
	output_charge.open("charge_lgad.dat", std::ios::out);
	output_time.open("time_lgad.dat", std::ios::out);
	output_gain.open("gain.dat", std::ios::out);
	output_signal.open("signal.txt", std::ios::out);
	output_y_real.open("yreal.txt");
	log.open("log_lgad", std::ios::out);
	log << "open" << "\n";
// signal simulation

	int numofpos = 0;
	double deltay = 2 * (gap + side) / numofpos;
	const unsigned int nEvents = 1;

	for (unsigned int k = 0; k < (numofpos + 1); ++k) {
	for (unsigned int j = 0; j < nEvents; ++j) {
		sensor.ClearSignal();
//		double y0 = 0.5 * (left + right) - numofpos * deltay + k * 2 * deltay;
//		double y0 = 0.5 * (left + right) + (RndmUniform() - 0.5) * (gap + side) * 4;
		double y0 = 0.5 * (left + right);
		output_y_real << y0 << "\n";
		std::cout << "y0 = " << y0 << "\n";
//		y0 += (RndmUniform() - 0.5) * (right - left) * 0.5; 
		double t0 = 0.1;
		track.NewTrack(0., y0, 0., t0, 0.1, 0., 0.);
		double xc = 0., yc = 0., zc = 0., tc = 0., ec = 0., dummy = 0.;
		int ne = 0, nh = 0;
		unsigned int nc = 0;
		unsigned int nesum = 0, nesum_real = 0, doublecheck = 0;
		std::cout << "position " << k << ", " << "No." << j << "\n";
		while (track.GetCluster(xc, yc, zc, tc, ne, nh, ec, dummy)) {
			++nc;
			nesum += ne; // primary electron number
			nesum_real += ne; // total electron number
			if (nc % 10 == 0) std::cout << "		Cluster " << nc << "\n";
			drift.DisablePlotting();
			if (plotDrift && RndmUniform() < 0.05) {
				drift.EnablePlotting(&vDrift);
			}
			double xe, ye, ze, te, ee, dxe, dye, dze;
			for (int i = 0; i < ne; ++i) {
				track.GetElectron(i, xe, ye, ze, te, ee, dxe, dye, dze);
//				drift.DriftElectron(xe, ye, ze, te); // for non-avalanche case
				drift.AvalancheElectronHole(xe, ye, ze, te);
				int size = drift.GetNumberOfElectronEndpoints();
				log << "size = " << size << "\n";
				nesum_real += size - 1;
				doublecheck += size; // another way of counting total electrons
				if (size > 1) log << "avalanche with " << size << "\n";
			}
		}

		std::cout << nesum << " electrons, " << nesum * ElementaryCharge << " fC.\n";
//		std::cout << nesum_real << " electrons, " << nesum_real * ElementaryCharge << " fC.\n";
		std::cout << doublecheck << " electrons, " << doublecheck * ElementaryCharge << " fC.\n";
		
		double gain = (double)nesum_real / nesum;
		output_gain << gain << "\n";

		output_charge << nesum_real * ElementaryCharge << "\n";
	  if (useTransfer)	sensor.ConvoluteSignals();

// write signals to files
		if (writeSignal) {
			for (unsigned int i = 0; i < nSignalBins; ++i) {
				const double t = (i + 0.5) * tStep;
				const double f_1 = sensor.GetSignal("pad1-delayed", i);
				const double f_2 = sensor.GetSignal("pad2-delayed", i);
				const double f_3 = sensor.GetSignal("pad3-delayed", i);
				const double f_4 = sensor.GetSignal("pad4-delayed", i);
				const double f_5 = sensor.GetSignal("pad5-delayed", i);
				output_signal << "  " << f_1 << "	" << f_2 << "	" << f_3 << "	" << f_4 << "	" << f_5 << "\n";
			}
		}

// write time info		
		if (writeTime) {
			std::cout << "write time" << "\n";
			double min = 999;
			double time = 0;
			double min_time = 0;
			// loop for peak value and corresponding time
			for (unsigned int i = 0; i < nSignalBins; ++i) {
				const double f = sensor.GetSignal(label1, i);
				const double t = (i + 0.5) * tStep;
				if (f < min) {
					min = f;
					min_time = t;
				}
			}
			// find the time stamp of half peak value in the rising period
			for (unsigned int i = 0; i < nSignalBins; ++i) {
				const double f = sensor.GetSignal(label6, i);
				const double t = (i + 0.5) * tStep;
				if (f > 0.5 * min && t < min_time) time = t;
			}
			output_time << time << "	" << min << "\n";
		}

		if (!plotSignal) continue;
		vSignal.PlotSignal("pad3-delayed", true, false, false, true, true);
		vSignal.PlotSignal("pad2-delayed", true, false, false, true, true);
		vSignal.PlotSignal("pad5-delayed", true, false, false, true, true);
		vSignal.PlotSignal("pad4-delayed", true, false, false, true, true);
		vSignal.PlotSignal("pad1-delayed", true, false, false, true, true);
	} // event loop over
	} // position loop over
	output_charge.close();
	output_time.close();
	output_gain.close();
	output_signal.close();
	output_y_real.close();
	log.close();

	if (plotDrift) {
		const bool twod = true;
		vDrift.Plot(twod, true);
	}

	app.Run(true); // if open, root will need to be closed manually

}
