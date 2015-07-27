#ifndef MARKER_H
#define MARKER_H

#include <vector>

#define LAST_CHROM	23

using namespace std;

class Marker {
  public:
    //////////////////////////////////////////////////////////////////
    // public static methods
    //////////////////////////////////////////////////////////////////

    static void readMarkerFile(char *markerFile);

    static int getNumMarkers() { return _allMarkers.size(); }
    static int getFirstMarkerNum(int chrom) { return _firstMarkerNum[chrom]; }
    static int getNumChromMarkers(int chrom)
		{ return _firstMarkerNum[chrom] - _firstMarkerNum[chrom - 1]; }

    static const Marker * getMarker(int num) { return _allMarkers[num]; }

    //////////////////////////////////////////////////////////////////
    // public methods
    //////////////////////////////////////////////////////////////////

    int getChrom() const { return _chrom; }
    double getMapPos() const { return _mapPos; }

  private:
    Marker(char *markerName, int chrom, double mapPos, int physPos,
	   char alleles[2]);


    // marker name (usually SNP rs id)
    char *_name;

    // chromosome
    int _chrom;

    // genetic position (map distance):
    double _mapPos;

    // physical position:
    int _physPos;

    // Major and minor alleles (respectively)
    char _alleles[2];

    //////////////////////////////////////////////////////////////////
    // private static variables
    //////////////////////////////////////////////////////////////////

    // List of all the markers read in:
    static vector<Marker *> _allMarkers;

    // Stores the number of marker loci on the corresponding chromosome number
    // 1..23
    static int _firstMarkerNum[LAST_CHROM + 1];
};


#endif // MARKER_H
