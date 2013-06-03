#ifndef __KML__Partitions__H
#define __KML__Partitions__H

#include <vector>
#include <iostream>
#include "Vector3D.h"

using namespace std; 

namespace KML{

class CParcel {

    friend class CPrivacyFineParcellatorHCP;
    friend class CFineParcellatorHCP;
    friend class CPartition;

    vector<float> timeSeries;
    vector<int> memberIds;
    int numVertices;
    int parcelLabel;
    RGBTYPE color;
    int centerVertexID;

};

class CPartition
{

    friend class CPrivacyFineParcellator;
    friend class CFineParcellator;

public:
    vector<float> principalVector;
    vector<int> parcelIds;
    vector<int> nbrPartitions;
    int numNbrParts; 
    float ratio;
    int numParcels;
    int partitionLabel;
    RGBTYPE color;
    VectorType sphereCoords;
    int centerVertexID;
    int numTimePoints;
    float area; 

public:

    CPartition();
    ~CPartition();

    CPartition& operator = (const CPartition& src);
    bool operator <(const CPartition& oth) const;
    bool operator >(const CPartition& oth) const;
    bool operator <= (const CPartition& oth) const;
    bool operator >= (const CPartition& oth) const;
    bool operator == (const CPartition& oth) const;
    bool operator!= (const CPartition& oth) const;
    void operator << (istream& strm);
    void operator >> (ostream& strm);

};

 ostream& operator << (ostream& out, CPartition& part);
 istream& operator >> (istream& in, CPartition& part); 

 
}; //end of namespace; 
#endif