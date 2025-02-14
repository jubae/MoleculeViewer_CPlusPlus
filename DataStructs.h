#ifndef DATASTRUCTS_H
#define DATASTRUCTS_H

#include <QDebug>
#include <QObject>
#include <QTOpenGL>

struct AtomPosition
{
    float x;
    float y;
    float z;
};

enum class EAtomColour
{
    White,
    Red,
    Green,
    Blue,
    Grey,
    Yellow,
    Orange,
    DarkRed,
    Black
};

enum class EAtomType
{
    Hydrogen = 0,
    Boron = 1,
    Carbon = 2,
    Nitrogen = 3,
    Oxygen = 4,
    Fluorine = 5,
    Silicon = 6,
    Phosphorus = 7,
    Sulphur = 8,
    Chlorine = 9,
    Copper = 10,
    Zinc = 11,
    Bromine = 12,
    Silver = 13,
    Iodine = 14,
    Gold = 15,
    Mercury = 16,
    Metal = 17,
    Unknown = 99
};

enum class EBondType
{
    Single = 1,
    Double = 2,
    Triple = 3,
    Unknown = 0
};

struct AtomData
{
    int AtomNumber;
    QVector3D AtomPosition;
    float AtomScale;
    QString AtomName;
    EAtomType AtomType;
};

struct BondData
{
    int BondNumber;
    int Atom1Number;
    int Atom2Number;
    EBondType BondType;
};

struct MoleculeData
{
    QVector<AtomData> Atoms;
    QVector<BondData> Bonds;
    QString MoleculeName;
};

struct AtomMeshData
{
    QVector<GLfloat> vertices;
    QVector<GLuint> indices;
    QVector<GLfloat> normals;
};


#endif //DATASTRUCTS_H