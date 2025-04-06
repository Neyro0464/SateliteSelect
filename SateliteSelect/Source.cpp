#include<iostream>
#include<string>
#include "Satselect.h"
#include <vector>

int main() {
	using namespace Sgpsdp;

	std::vector<NORAD_DATA> Satellites;
	const std::string file = "TLE.txt";
	Station st;

	st.geo.Lat = 56.14717;
	st.geo.Lon = 37.76275;
	st.geo.Alt = 100;

	st.lim.minAzm = 120.0;
	st.lim.maxAzm = 200.0;
	st.lim.minElv = 10.0;
	st.lim.maxElv = 80.0;
	st.lim.timeMinObserveSec = 20;

	SetelliteSelect tmp(file, st);
	Satellites = tmp.GetSatArray();
	tmp.showSat(Satellites);
	return 0;
}