#ifndef SATSELECT_H
#define SATSELECT_H

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <sstream>
#include "Sgpsdp.h"
#include <cstring>

namespace Sgpsdp {

	struct CoordDecart {
		double x = 0;
		double y = 0;
		double z = 0;
	};
	struct CoordTopocentric {
		double azm = 0;
		double elv = 0;
	};
	struct CoordGeodetic {
		double Lat = 0;
		double Lon = 0;
		double Alt = 0;
	};
	struct StationVision {

		double minAzm = 0;
		double maxAzm = 0;
		double minElv = 0;
		double maxElv = 0;
		double timeMinObserveSec = 0;
	};
	struct Station {
		CoordDecart dec{};
		CoordGeodetic geo{};
		StationVision lim{};
	};
	struct COORDS {
		CoordTopocentric topo{};
		CoordGeodetic geo{};
	};
	struct NORAD_DATA {
		std::string name{};
		unsigned int noradNumber = 0;
		std::string onTime{};
		COORDS coords{};
		double time = 0;
		bool dirPositive{};

	};
	class SetelliteSelect {
	private:
		std::vector<CSGP4_SDP4> Satellites;
		Station station;
		std::string filename;

		double DegToRad(double arg);
		double RadToDeg(double arg);
		CoordDecart ConvertGEOtoDecart(CoordGeodetic& Object);
		CoordDecart Transferring(CoordDecart object, CoordGeodetic Station);

		bool VelocityCheck(Station station, CSGP4_SDP4 SatelliteModel);
		bool IntersectCheck(CoordDecart& Satellite_1, Station station, double& theta, double& fi);
		bool getDiraction(Station station, CSGP4_SDP4& SatelliteModel);
		double TimeIntersect(CSGP4_SDP4& SatelliteModel, Station& station, double& theta, double& fi);
		CoordDecart SatellitePos(CSGP4_SDP4 SGP, CoordGeodetic& SatelliteLLA_1, CoordDecart Satellite_1, double time);
		void dataPrepare();


	public:

		SetelliteSelect(const std::string& filename, Station& param);
		~SetelliteSelect() = default;

		// Enter from console station coordinates (not in use)
		bool SetStationPos(Station& station);
		// Enter from console station params (not in use)
		bool SetFilter(Station& station);

		/*
		Function return final satellite array
		*/
		const std::vector<NORAD_DATA> GetSatArray();
		/*
		Function output final satellite array to the console
		*/
		void showSat(const std::vector<NORAD_DATA>& result);


	};

}

#endif // SGP4_H
