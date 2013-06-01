#include "kmlPartitions.h"
#include "iostream"
#include "fstream"

using namespace std; 


namespace KML{


CPartition::CPartition()
{
    this->centerVertexID=-1;
    this->color= RGBTYPE ( 0,0,0 );
    this->numParcels=0;
    this->numTimePoints=0;
    this->ratio=1;
    this->partitionLabel=-1;
    this->sphereCoords=VectorType ( 0,0,0 );
    this->area=0;
}
CPartition::~CPartition()
{

}

bool CPartition::operator!= ( const KML::CPartition& oth ) const
{
    if ( partitionLabel !=oth.centerVertexID || centerVertexID != oth.centerVertexID || color!= oth.color
            || sphereCoords != oth.sphereCoords || numTimePoints != oth.numTimePoints || area != oth.area ) {
        return false;
    }

    for ( int idx = 0; idx < numTimePoints ; ++idx ) {
        if ( principalVector[idx] != oth.principalVector[idx] ) {
            return false;
        }
    } // end for loop::idx

    return true;
}
void CPartition::operator<< ( istream& strm )
{
    strm.read ( ( char* ) &partitionLabel, sizeof ( int ) );
    strm.read ( ( char* ) &centerVertexID, sizeof ( int ) );
    strm.read ( ( char* ) & ( color.x ), 3*sizeof ( float ) );
    strm.read ( ( char* ) & ( sphereCoords.x ), 3*sizeof ( float ) );
    strm.read ( ( char* ) &numTimePoints, sizeof ( int ) );
    principalVector.resize ( numTimePoints,0 );
    strm.read ( ( char* ) &principalVector[0],  numTimePoints*sizeof ( float ) );
    strm.read ( ( char* ) &area,  sizeof ( float ) );
    strm.read ( ( char* ) &ratio,  sizeof ( float ) );
    strm.read ( ( char* ) &numNbrParts,sizeof ( int ) );
    nbrPartitions.resize ( numNbrParts,0 );
    strm.read ( ( char* ) & ( nbrPartitions[0] ),numNbrParts*sizeof ( int ) );
    strm.read ( ( char* ) &numParcels,sizeof(int));
    parcelIds.resize(numParcels,0);    
    strm.read ( ( char* ) & ( parcelIds[0] ),numParcels*sizeof ( int ) );

}
void CPartition::operator>> ( ostream& strm )
{
    strm.write ( ( char* ) &partitionLabel, sizeof ( int ) );
    strm.write ( ( char* ) &centerVertexID, sizeof ( int ) );
    strm.write ( ( char* ) & ( color.x ), 3*sizeof ( float ) );
    strm.write ( ( char* ) & ( sphereCoords.x ), 3*sizeof ( float ) );
    strm.write ( ( char* ) &numTimePoints, sizeof ( int ) );
    strm.write ( ( char* ) & ( principalVector[0] ), numTimePoints*sizeof ( float ) );
    strm.write ( ( char* ) & ( area ), sizeof ( float ) );
    strm.write ( ( char* ) & ( ratio ), sizeof ( float ) );
    strm.write ( ( char* ) &numNbrParts,sizeof ( int ) );
    strm.write ( ( char* ) & ( nbrPartitions[0] ),numNbrParts*sizeof ( int ) );
    strm.write ( ( char* ) &numParcels,sizeof(int));
    strm.write ( ( char* ) & ( parcelIds[0] ),numParcels*sizeof ( int ) );
}

CPartition& CPartition::operator= ( const KML::CPartition& src )
{
    if ( this == &src ) {
        return *this;
    } else {
        this->centerVertexID= src.centerVertexID;
        this->color=src.color;
        this->numParcels=src.numParcels;
        this->parcelIds=src.parcelIds;
        this->partitionLabel=src.partitionLabel;
        this->principalVector=src.principalVector;
        this->ratio=src.ratio;
        this->sphereCoords=src.sphereCoords;
    }
}
bool CPartition::operator== ( const KML::CPartition& oth ) const
{
    if ( this->numParcels == oth.numParcels && this->ratio == oth.ratio ) {
        return true;
    } else {
        return false;
    }
}

bool CPartition::operator< ( const KML::CPartition& oth ) const
{

    float v = float ( numParcels ) * ( 1.f-ratio );
    float voth = float ( oth.numParcels ) * ( 1-oth.ratio );
    if ( v > voth ) {
        return true;
    }
    return false;
    if ( ratio < oth.ratio ) {
        return true;
    }
    return false;
}

bool CPartition::operator<= ( const KML::CPartition& oth ) const
{
    return this->operator== ( oth ) || this->operator< ( oth );
}
bool CPartition::operator> ( const KML::CPartition& oth ) const
{
    return !this->operator<= ( oth )  ;
}
bool CPartition::operator>= ( const KML::CPartition& oth ) const
{
    return !this->operator< ( oth );
}

ostream& operator << ( ostream& out, CPartition& part )
{
    part.operator>> ( out );
    return out;

}
istream& operator >> ( istream& in, CPartition& part )
{
    part.operator<< ( in );
    return in;
}

};//end of namespace KML: 
