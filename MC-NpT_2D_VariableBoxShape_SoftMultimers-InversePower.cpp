#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "MTGenerator.h"

/*State of optional functions:  (remember to always update text below after changing this file)
 -adjustAngles - ACTIVE
 -adjustOddRowsTranslation - INACTIVE
 -rotationsAnalysis - ACTIVE
 -Rectangular volume moves (not romboidal) - INACTIVE
*/

int N,gaps,activeN,loadedConfiguration,loadType=0,loadedSetStartGenerator,loadedSetGenerator,iterationsNumber,
    growing,multimerN,countCollidingPairs,ODFLength,OCFMode,skipFirstIteration,saveConfigurations,
    useSpecificDirectory,useFileToIterate,fileIterateIterationsNumber=0,actIteration=0,multiplyArgument,
    onlyMath[2]={0,0},initMode,autoEqLength,neighUpdatingFrequency;
long cyclesOfEquilibration,cyclesOfMeasurement,timeEq=0,timeMe=0,timeMath=0,intervalSampling,intervalOutput,intervalResults,intervalOrientations,savedConfigurationsInt,generatorStartPoint=0;
double maxDeltaR,desiredAcceptanceRatioR,desiredAcceptanceRatioV,
       startMinPacFrac,startMaxPacFrac,minArg,maxArg,loadedArg,
       intervalMin[10],intervalMax[10],intervalDelta[10],
       startArg[2],deltaR,deltaPhi,deltaV=0.1, //*d^**\sigma
       TReduced,invPotPower,multimerS,multimerD,randomStartStep[2],detBoxMatrix,
       neighRadius,neighRadius2,multiplyFactor,pressureRealOfNotFluid,
       iterationTable[1000][3],pi=M_PI; //Problem typu-Pi (program przestawał działać dla wyższej precyzji Pi) to problem precyzji double. AKURAT dla 0.500001 się on zaczyna (przy 0.500010 jeszcze wszystko działa dobrze). Chodzi o to, że dla programu w C++ d=0.500001 to tak naprawdę: 0.500001000000000023 (np.), a dla takiej średnicy, procedury zwracają już INNE absoluteMinimum i minDistance. Normalnie cyfry tak 'dalekie' mogą być pominięte, ale przy d=0.500001 zaczyna już to odgrywać rolę. Zmniejszenie precyzji PI tak naprawdę ZMNIEJSZAŁO PI (ucięcie ostatnich liczb rozwinięcia) co jakoś przekładało się na to, że program zwracał 'wyglądające na dobre' wyniki.
double L,C,ROkreguOpisanego,hTrojkataWielokata,absoluteMinimum,absoluteMinimum2,minDistance,maxDistance,maxDistanceCore,VcpPerParticle,dr[3];
char buffer[200]="",bufferN[20],bufferGaps[20],bufferG[5],bufferMN[20],bufferMS[20],bufferMD[20],bufferFolderIndex[5],bufferTReduced[20],
     resultsFileName[200]="ResultsSummary.txt",
     excelResultsFileName[200]="ExcelResultsSummary.txt",
     configurationsFileName[200]="Configurations",
     orientationsFileName[200]="Orientations",
     orientatCorrelFunFileName[200]="OrientatCorrelFun",
     orientationsResultsFileName[200]="OrientatRes",
     configurationsListFileName[200]="ConfigurationsList.txt",
     loadConfigurationsFileName[200]="Configurations",
     loadedJOBID[50]="j-none";
/////////////////  PARTICLE functions{
typedef struct {
    double r[2], normR[2], atomsR[20][2];  //x,y; atomsR - względne pozycje atomów, dla max multimerN=20
    double phi;   //kąt mierzony od kierunku x
    int neighbours[50], neighCounter;   //uwaga na przepełnienie neighbours[] - wówczas przy tworzeniu listy sąsiadów program może 'pisać' w innych zmiennych, np. w r[] zamiast informować (już tak bywało)
    double neighEnergyFactor[50];
} particle;


int getIndexFromList (int *list, int listLength, int searchedElement) {
    int indexInList=-1;
    for (int i=0;i<listLength;i++) if (list[i]==searchedElement) {
        indexInList=i;
        break;
    }
    return indexInList;
}

void getParticlesDistanceSquared (particle *p1, particle *p2, double boxMatrix[2][2]) {
    double normalizedRX=p1->normR[0]-p2->normR[0],
           normalizedRY=p1->normR[1]-p2->normR[1],
           rx=p1->r[0]-p2->r[0],
           ry=p1->r[1]-p2->r[1];
    rx-=round(normalizedRX)*boxMatrix[0][0]+round(normalizedRY)*boxMatrix[0][1];
    ry-=round(normalizedRX)*boxMatrix[1][0]+round(normalizedRY)*boxMatrix[1][1];
    dr[0]=rx; dr[1]=ry; dr[2]=rx*rx+ry*ry;
}

void computeAtomsPositions (particle *p) {
    for (int i=0;i<multimerN;i++) {
        double discAngle=p->phi+i*2*C;
        p->atomsR[i][0]=cos(discAngle)*ROkreguOpisanego;
        p->atomsR[i][1]=sin(discAngle)*ROkreguOpisanego;
    }
}

/*void adjustNeighRadiusMod (double *factor, int mod) {
    neighRadiusMod+=*factor*mod;
    if (neighRadiusMod>0) {
        if (neighRadiusMod>=2) *factor=0.000000000000001;
        neighRadius=neighRadiusMod*maxDistance;
        neighRadius2=neighRadius*neighRadius;
    } else {
        neighRadiusMod-=*factor*mod;
        *factor*=0.5;
    }
}

void updateNeighbourList (particle *particles, double boxMatrix[2][2], int numberOfClosestNeighbours) {
    int neighRadiusChange[2]={0,0}; double neighRadiusChangeFactor=0.1,oldNeighRadiusMod=neighRadiusMod; bool exception;
    for (int i=0;i<activeN;i++) particles[i].neighCounter=0;
    for (int i=0;i<activeN-1;i++) do {exception=false; int bufferNumberOfClosestNeighbours=numberOfClosestNeighbours;
        int neighCounterBuffer=particles[i].neighCounter;
        if (neighCounterBuffer>=numberOfClosestNeighbours) bufferNumberOfClosestNeighbours=neighCounterBuffer+1; //kazda czastka ma miec prawo znalezienia PRZYNAJMNIEJ jednej najblizszej SPOSROD POZOSTALYCH (nizej w tabeli)
        for (int j=i+1;j<activeN;j++) {
            getParticlesDistanceSquared(particles[i],particles[j],boxMatrix);
            if (dr[2]<neighRadius2) particles[i].neighbours[particles[i].neighCounter++]=j;
            if (particles[i].neighCounter>bufferNumberOfClosestNeighbours) break;
        }
        if (particles[i].neighCounter!=bufferNumberOfClosestNeighbours && neighRadiusChangeFactor>=0.00000000000001) {
            neighRadiusChange[1]=neighRadiusChange[0];
            neighRadiusChange[0]=particles[i].neighCounter>bufferNumberOfClosestNeighbours?-1:1;
            if (neighRadiusChange[0]*neighRadiusChange[1]<0) neighRadiusChangeFactor*=0.5;
            adjustNeighRadiusMod(&neighRadiusChangeFactor,neighRadiusChange[0]);
            particles[i].neighCounter=neighCounterBuffer;
        } else {
            for (int j=neighCounterBuffer;j<particles[i].neighCounter;j++)
                particles[particles[i].neighbours[j]].neighbours[particles[particles[i].neighbours[j]].neighCounter++]=i;
            if (particles[i].neighCounter!=numberOfClosestNeighbours) {
                printf("Couldn't get desired number of closest neighbours. Got: %d for particle: %d\n",particles[i].neighCounter,i);
                exception=true;
            }
            neighRadiusChange[0]=0; neighRadiusChangeFactor=0.1;
        }
    } while (particles[i].neighCounter!=numberOfClosestNeighbours && !exception);
    if (oldNeighRadiusMod!=neighRadiusMod) {
        printf("Adjusted neighRadiusMod: %.17E.\n",neighRadiusMod); fflush(stdout);
    }
}*/

void adjustNeighRadius (double volume) {
    neighRadius=1.6*sqrt(volume)/sqrt(N);
    neighRadius2=neighRadius*neighRadius;
}

void updateNeighbourList (particle *particles, double boxMatrix[2][2], double volume) {
    adjustNeighRadius(volume);
    for (int i=0;i<activeN;i++) particles[i].neighCounter=0;
    for (int i=0;i<activeN-1;i++) for (int j=i+1;j<activeN;j++) {
        getParticlesDistanceSquared(&particles[i],&particles[j],boxMatrix);
        if (dr[2]<neighRadius2) {
            particles[i].neighbours[particles[i].neighCounter++]=j;
            particles[j].neighbours[particles[j].neighCounter++]=i;
        }
    }
}

double normalizeAngle (double phi) {//inaczej niż w Mathematice (loopedAngle[phi_]:=Mod[phi+C,2*C]-C), bo Mod[a,b] działa inaczej. Dla a>0 i b>0 nie ma rożnicy, ale dla a<0 i b>0 JEST. Mathematica zwraca wartości dodatnie (Mod[-1,10]=9), a C++ nie (fmod(-1,10)=-1).
    phi=fmod(phi+C,2*C);
    return phi<0?phi+C:phi-C;
}

void checkSinglePeriodicBoundaryConditions (particle *p, double boxMatrix[2][2]) {
    for (int j=0;j<2;j++) {
        for (int i=0;i<2;i++) p->r[i]-=floor(p->normR[j])*boxMatrix[i][j];
        p->normR[j]=fmod(p->normR[j],1); if (p->normR[j]<0) p->normR[j]++;
    }
    //particle->phi=normalizeAngle(particle->phi);    //rotationsAnalysis(comment) 1/1
}

void checkPeriodicBoundaryConditions (particle *particles, double boxMatrix[2][2]) {
    for (int i=0;i<activeN;i++)
        checkSinglePeriodicBoundaryConditions(&particles[i],boxMatrix);
}

double minimalDistanceAnalyticalMethodHCM (int indeksPowierzchni, double aAngle, double bAngle, double a[4], double b[4]) {
    double angleA=aAngle+a[indeksPowierzchni],angleB=bAngle+b[indeksPowierzchni],
           buffer=sin(angleA)+sin(angleB);
    return ROkreguOpisanego*(cos(angleA)+cos(angleB))+sqrt(multimerD*multimerD-ROkreguOpisanego*ROkreguOpisanego*buffer*buffer);
}

double getMinimalDistanceAnalyticalMethodForEvenHCM (double aAngle, double bAngle) {
    double a[4]={C,-C,C,-C}, b[4]={-C,-C,C,C};
    if (aAngle<-absoluteMinimum) {
        double ARCSIN=asin(sin(aAngle+C)*L/(L*cos(C)+sqrt(4.0-L*L)*sin(C))),
               intersectionCurve13=-ARCSIN;
        if (bAngle>intersectionCurve13) return minimalDistanceAnalyticalMethodHCM(0,aAngle,bAngle,a,b);
        else return minimalDistanceAnalyticalMethodHCM(2,aAngle,bAngle,a,b);
    } else if (aAngle<0) {
        double ARCSIN=asin(sin(aAngle)/L*(L*cos(C)+sqrt(4.0-L*L)*sin(C))),
               intersectionCurve34=-C-ARCSIN;
        double intersectionCurve14=aAngle;
        if (bAngle>intersectionCurve14) return minimalDistanceAnalyticalMethodHCM(0,aAngle,bAngle,a,b);
        else if (bAngle>intersectionCurve34) return minimalDistanceAnalyticalMethodHCM(3,aAngle,bAngle,a,b);
        else return minimalDistanceAnalyticalMethodHCM(2,aAngle,bAngle,a,b);
    } else if (aAngle<absoluteMinimum) {
        double ARCSIN=asin(sin(aAngle)/L*(L*cos(C)+sqrt(4.0-L*L)*sin(C))),
               intersectionCurve12=C-ARCSIN;
        double intersectionCurve14=aAngle;
        if (bAngle>intersectionCurve12) return minimalDistanceAnalyticalMethodHCM(1,aAngle,bAngle,a,b);
        else if (bAngle>intersectionCurve14) return minimalDistanceAnalyticalMethodHCM(0,aAngle,bAngle,a,b);
        else return minimalDistanceAnalyticalMethodHCM(3,aAngle,bAngle,a,b);
    } else {
        double ARCSIN=asin(sin(aAngle-C)*L/(L*cos(C)+sqrt(4.0-L*L)*sin(C))),
               intersectionCurve24=-ARCSIN;
        if (bAngle>intersectionCurve24) return minimalDistanceAnalyticalMethodHCM(1,aAngle,bAngle,a,b);
        else return minimalDistanceAnalyticalMethodHCM(3,aAngle,bAngle,a,b);
    }
}

double getMinimalDistanceAnalyticalMethodForOddHCM (double aAngle, double bAngle) {
    double a[4]={C,-C,C,-C}, b[4]={-2*C,0,0,2*C};
    if (aAngle<-absoluteMinimum) {
        double ARCSIN = asin(sin(aAngle+C)*L/(L*cos(C)+sqrt(4.0-L*L)*sin(C))),
               intersectionCurve13 = C-ARCSIN;
        if (bAngle>intersectionCurve13) return minimalDistanceAnalyticalMethodHCM(0,aAngle,bAngle,a,b);
        else return minimalDistanceAnalyticalMethodHCM(2,aAngle,bAngle,a,b);
    } else if (aAngle<0) {
        double ARCSIN = asin(sin(aAngle)/L*(L*cos(C)+sqrt(4.0-L*L)*sin(C))),
               intersectionCurve23 = -ARCSIN;
        double intersectionCurve12 = aAngle+C;
        if (bAngle>intersectionCurve12) return minimalDistanceAnalyticalMethodHCM(0,aAngle,bAngle,a,b);
        else if (bAngle>intersectionCurve23) return minimalDistanceAnalyticalMethodHCM(1,aAngle,bAngle,a,b);
        else return minimalDistanceAnalyticalMethodHCM(2,aAngle,bAngle,a,b);
    } else if (aAngle<absoluteMinimum) {
        double ARCSIN = asin(sin(aAngle)/L*(L*cos(C)+sqrt(4.0-L*L)*sin(C))),
               intersectionCurve23 = -ARCSIN;
        double intersectionCurve34 = aAngle-C;
        if (bAngle>intersectionCurve23) return minimalDistanceAnalyticalMethodHCM(1,aAngle,bAngle,a,b);
        else if (bAngle>intersectionCurve34) return minimalDistanceAnalyticalMethodHCM(2,aAngle,bAngle,a,b);
        else return minimalDistanceAnalyticalMethodHCM(3,aAngle,bAngle,a,b);
    } else {
        double ARCSIN = asin(sin(aAngle-C)*L/(L*cos(C)+sqrt(4.0-L*L)*sin(C))),
               intersectionCurve24 = -C-ARCSIN;
        if (bAngle>intersectionCurve24) return minimalDistanceAnalyticalMethodHCM(1,aAngle,bAngle,a,b);
        else return minimalDistanceAnalyticalMethodHCM(3,aAngle,bAngle,a,b);
    }
}

int checkOverlapsAnalyticalMethodHCM (double dr, double aAngle, double bAngle) {
    int overlap=0;

    aAngle=normalizeAngle(aAngle+C);
    bAngle=normalizeAngle(bAngle+C);

    if (multimerN%2==0) {
        if (dr<getMinimalDistanceAnalyticalMethodForEvenHCM(aAngle,bAngle)) overlap=1;
    } else if (dr<getMinimalDistanceAnalyticalMethodForOddHCM(aAngle,bAngle)) overlap=1;
    return overlap;
}

int checkMoleculeInteriorsOverlap (double dr, double aAngle, double bAngle) { //energy<0 (-1) -> overlap
    int energy=0;

    aAngle=normalizeAngle(aAngle+C);
    bAngle=normalizeAngle(bAngle+C);
    double absAAngle=fabs(aAngle),absBAngle=fabs(bAngle);

    if (absBAngle>absAAngle) {
        double mod=C-absBAngle;
        if(dr<hTrojkataWielokata/cos(aAngle)+ROkreguOpisanego*sin(mod)*tan((bAngle>0)?-aAngle:aAngle)+ROkreguOpisanego*cos(mod)) energy=-1;
    } else {
        double mod=C-absAAngle;
        if(dr<hTrojkataWielokata/cos(bAngle)+ROkreguOpisanego*sin(mod)*tan((aAngle>0)?-bAngle:bAngle)+ROkreguOpisanego*cos(mod)) energy=-1;
    }
    return energy;
}

int createRandomGaps (particle *particles, double boxMatrix[2][2], double volume) { //TODO: 1) głupie swapowanie - po co przesuwać cząstki w liście, jak wystarczy po prostu 'przerzucać' gapy na (aktualny) koniec listy i zmniejszać jej wielkość (ucinając automatycznie ten ostatni element); 2) nie zamienia wszystkich pól obiektów cząstek (przestarzałe);
    printf("Creating %d random gaps... ",gaps); fflush(stdout);
    adjustNeighRadius(volume);
    int gapsIndexes[gaps], attempt=0, innerAttempt=0, allReady=0;
    do {
        attempt++;
        for (int i=0;i<gaps;i++) {
            gapsIndexes[i] = (int)(MTRandom0to1(randomStartStep)*N);
            for (int j=0;j<i;j++) {
                getParticlesDistanceSquared(&particles[gapsIndexes[i]],&particles[gapsIndexes[j]],boxMatrix);
                if (dr[2]<neighRadius2) {
                    i--;
                    innerAttempt++;
                    break;
                }
                if (j+1==i) innerAttempt=0;
            }
            if (innerAttempt>100000) break;
            if (innerAttempt==0 && i+1==gaps) allReady=1;
        }
    } while (!allReady && attempt<10000000);

    if (attempt>=10000000) {
        printf("ERROR: Couldn't create %d random gaps in %d steps.\n",gaps,attempt);
        return 0;
    } else {
        bool change; do {
            change=false;
            for (int i=gaps-1;i>0;i--) {
                if (gapsIndexes[i-1]>gapsIndexes[i]) {
                    int buffer=gapsIndexes[i];
                    gapsIndexes[i]=gapsIndexes[i-1]; gapsIndexes[i-1]=buffer;
                    change=true;
                }
            }
        } while (change);
        int actualGapIndex=0;
        for (int i=0;i<activeN;i++) {
            if (actualGapIndex<gaps) while (i+actualGapIndex==gapsIndexes[actualGapIndex]) {
                actualGapIndex++;
                if (actualGapIndex==gaps) break;
            }
            for (int j=0;j<2;j++) {
                particles[i].r[j]=particles[i+actualGapIndex].r[j];
                particles[i].normR[j]=particles[i+actualGapIndex].normR[j];
            }
        }
        printf("done\n"); fflush(stdout);
        return 1;
    }
}

int initPositions (particle *particles, double boxMatrix[2][2], double matrixOfParticlesSize[2], int n[2], double matrixCellXY[6][6][2], double matrixCellPhi[6][6], double pacFrac, double volume) {
    if (generatorStartPoint==0) {
        printf("Setting start position of p-random number generator to actual CPU time (for INIT PROCEDURES)...\n");
        InitRandomMT();
    } else {
        printf("Setting start position of p-random number generator to %ld (for INIT PROCEDURES)...\n",generatorStartPoint);
        InitMT((unsigned int)generatorStartPoint);
    }

    double mod=sqrt(N/matrixOfParticlesSize[0]/matrixOfParticlesSize[1]), interval[2][2], actualPosition[2];
    for (int i=0;i<2;i++) for (int j=0;j<2;j++) interval[i][j]=boxMatrix[i][j]/matrixOfParticlesSize[j]/mod*n[j];
    int rowCounter=0, columnCounter=0;
    for (int i=0;i<N;i++) {
        int cellNumber[2][2]={{columnCounter/n[0],rowCounter/n[1]},{columnCounter%n[0],rowCounter%n[1]}}; //cellNumber[0/1][X/Y]: 0-numer komorki, 1-kolumna/rzad W komorce
        actualPosition[0]=cellNumber[0][0]*interval[0][0]+cellNumber[0][1]*interval[0][1]+matrixCellXY[cellNumber[1][0]][cellNumber[1][1]][0]*sqrt(pacFrac);
        actualPosition[1]=cellNumber[0][1]*interval[1][1]+cellNumber[0][0]*interval[1][0]+matrixCellXY[cellNumber[1][0]][cellNumber[1][1]][1]*sqrt(pacFrac);

        for (int j=0;j<2;j++) particles[i].r[j]=actualPosition[j];
        particles[i].normR[0]=(boxMatrix[1][1]*particles[i].r[0]-boxMatrix[0][1]*particles[i].r[1])/detBoxMatrix;
        particles[i].normR[1]=-(boxMatrix[1][0]*particles[i].r[0]-boxMatrix[0][0]*particles[i].r[1])/detBoxMatrix;
        particles[i].phi=matrixCellPhi[cellNumber[1][0]][cellNumber[1][1]];
        computeAtomsPositions(&particles[i]);

        columnCounter++;
        if (columnCounter*1.000001>=matrixOfParticlesSize[0]*mod) {
            rowCounter++;
            columnCounter=0;
        }
    }
    if (gaps>0) return createRandomGaps(particles,boxMatrix,volume);
    else return 1;
}

double getEnergy (particle *p1, particle *p2) {
    double energy=0;
    for (int i=0;i<multimerN;i++) for (int j=0;j<multimerN;j++) {
        double xDistance=dr[0]+p1->atomsR[i][0]-p2->atomsR[j][0], yDistance=dr[1]+p1->atomsR[i][1]-p2->atomsR[j][1],
               dR=sqrt(xDistance*xDistance+yDistance*yDistance);
        energy+=pow(multimerD/dR,invPotPower);
    }
    return energy;
}

double getEnergyAll (particle *particles, double boxMatrix[2][2], int checkType) {  //energy<0 (-1) -> core overlap, checkType: 0-both, 1-only mol overlap, 2-only energy calculation
    double energy=0;
    if (checkType!=2) for (int i=0;i<activeN-1;i++) for (int j=i+1;j<activeN;j++) {
        getParticlesDistanceSquared(&particles[i],&particles[j],boxMatrix);
        double drRoot=sqrt(dr[2]);
        if (drRoot<maxDistanceCore) {
            double gamma=atan(dr[1]/dr[0]),
                   aAngle=particles[j].phi-gamma,
                   bAngle=particles[i].phi-gamma;
            if (multimerN%2!=0) {
                //rozważanie, która molekuła jest 'po lewej', a która 'po prawej'
                if (dr[0]>0) bAngle-=C;
                else aAngle-=C;
            }
            energy=checkMoleculeInteriorsOverlap(drRoot,aAngle,bAngle);
            if (energy<0) {
                i=activeN; j=activeN; break;
            }
        }
    }
    if (checkType!=1 && energy>=0) for (int i=0;i<activeN;i++) for (int j=0;j<particles[i].neighCounter;j++) {
        if (i<particles[i].neighbours[j]) {
            getParticlesDistanceSquared(&particles[i],&particles[particles[i].neighbours[j]],boxMatrix);
            particles[i].neighEnergyFactor[j]=getEnergy(&particles[i],&particles[particles[i].neighbours[j]]);
            particles[particles[i].neighbours[j]].neighEnergyFactor[getIndexFromList(particles[particles[i].neighbours[j]].neighbours,particles[particles[i].neighbours[j]].neighCounter,i)]=particles[i].neighEnergyFactor[j];
            energy+=particles[i].neighEnergyFactor[j];
        }
    }
    return energy;
}

void adjustAngles (particle *particles, double boxMatrix[2][2]) {
    int tryNumber=-1,energy;
    printf("Angle adjusting... "); fflush(stdout);
    do {
        tryNumber++;
        energy=(int)getEnergyAll(particles,boxMatrix,1);
        if (energy<0) for (int k=0;k<activeN;k++) {particles[k].phi+=0.00001; computeAtomsPositions(&particles[k]);}
    } while (energy<0);
    printf("Adjusted after: %d approach.\n",tryNumber); fflush(stdout);
}

void adjustOddRowsTranslation (particle *particles, double boxMatrix[2][2], int partInRow) {
    int tryNumber=-1,energy;
    printf("Odd rows translation adjusting... "); fflush(stdout);
    do {
        tryNumber++;
        energy=(int)getEnergyAll(particles,boxMatrix,1);
        if (energy<0) for (int k=0;k<activeN;k++) if ((k/partInRow)%2==1) {
            particles[k].r[0]+=0.00000001*multimerD;  //jednostka odległości to d^**\sigma
            particles[k].normR[0]=(boxMatrix[1][1]*particles[k].r[0]-boxMatrix[0][1]*particles[k].r[1])/detBoxMatrix;
            particles[k].normR[1]=-(boxMatrix[1][0]*particles[k].r[0]-boxMatrix[0][0]*particles[k].r[1])/detBoxMatrix;
            checkSinglePeriodicBoundaryConditions(&particles[k],boxMatrix);
        }
    } while (energy<0);
    printf("Adjusted after: %d approach.\n",tryNumber); fflush(stdout);
}

double getVicinityEnergyChange (int *result, particle *particles, particle *dispPart, int index, double boxMatrix[2][2]) {
    double deltaEnergy=0;
    for (int i=0;i<particles[index].neighCounter;i++) {
        getParticlesDistanceSquared(&particles[particles[index].neighbours[i]],dispPart,boxMatrix);
        double drRoot=sqrt(dr[2]);
        if (drRoot<maxDistanceCore) {
            double gamma=atan(dr[1]/dr[0]),
                   aAngle=dispPart->phi-gamma,
                   bAngle=particles[particles[index].neighbours[i]].phi-gamma;
            if (multimerN%2!=0) {
                //rozważanie, która molekuła jest 'po lewej', a która 'po prawej'
                if (dr[0]>0) bAngle-=C;
                else aAngle-=C;
            }
            if (checkMoleculeInteriorsOverlap(drRoot,aAngle,bAngle)<0) {
                i=particles[index].neighCounter; *result=0; break;
            }
        }
    }
    if (*result==1) for (int i=0;i<particles[index].neighCounter;i++) {
        getParticlesDistanceSquared(&particles[particles[index].neighbours[i]],dispPart,boxMatrix);
        dispPart->neighEnergyFactor[i]=getEnergy(&particles[particles[index].neighbours[i]],dispPart);
        deltaEnergy+=dispPart->neighEnergyFactor[i]-particles[index].neighEnergyFactor[i];
    }
    return deltaEnergy;
}

int attemptToDisplaceAParticle (particle *particles, int index, double boxMatrix[2][2], double *totalInteractionEnergy) {
    int result=1;
    particle displacedParticle;
    for (int i=0;i<2;i++) displacedParticle.r[i]=particles[index].r[i]+(MTRandom0to1(randomStartStep)-0.5)*deltaR;
    displacedParticle.normR[0]=(boxMatrix[1][1]*displacedParticle.r[0]-boxMatrix[0][1]*displacedParticle.r[1])/detBoxMatrix;
    displacedParticle.normR[1]=-(boxMatrix[1][0]*displacedParticle.r[0]-boxMatrix[0][0]*displacedParticle.r[1])/detBoxMatrix;
    displacedParticle.phi=particles[index].phi+(MTRandom0to1(randomStartStep)-0.5)*deltaPhi;
    computeAtomsPositions(&displacedParticle);
    double deltaEnergy=getVicinityEnergyChange(&result,particles,&displacedParticle,index,boxMatrix);

    if (result==1) {
        double arg=-deltaEnergy/TReduced; //T*=kT/eps [eps - jednostka energii (w domyśle =1)]
        if (MTRandom0to1(randomStartStep)>exp(arg)) result=0;
        else {
            *totalInteractionEnergy+=deltaEnergy;
            for (int i=0;i<2;i++) {
                particles[index].r[i]=displacedParticle.r[i];
                particles[index].normR[i]=displacedParticle.normR[i];
                for (int j=0;j<multimerN;j++) particles[index].atomsR[j][i]=displacedParticle.atomsR[j][i];
            }
            particles[index].phi=displacedParticle.phi;
            for (int i=0;i<particles[index].neighCounter;i++) {
                particles[index].neighEnergyFactor[i]=displacedParticle.neighEnergyFactor[i];
                particles[particles[index].neighbours[i]].neighEnergyFactor[getIndexFromList(particles[particles[index].neighbours[i]].neighbours,particles[particles[index].neighbours[i]].neighCounter,index)]=displacedParticle.neighEnergyFactor[i];
            }
            checkSinglePeriodicBoundaryConditions(&particles[index],boxMatrix);
        }
    }
    return result;
}

void cloneParticlesForSpecificBoxMatrix (particle* clonedParticles, particle *particles, double boxMatrix[2][2]) {
    for (int i=0;i<activeN;i++) {
        for (int j=0;j<2;j++) {
            clonedParticles[i].normR[j]=particles[i].normR[j];
            for (int k=0;k<multimerN;k++) clonedParticles[i].atomsR[k][j]=particles[i].atomsR[k][j];
        }
        for (int j=0;j<2;j++) clonedParticles[i].r[j]=boxMatrix[j][0]*particles[i].normR[0]+boxMatrix[j][1]*particles[i].normR[1];
    }
}

int attemptToChangeVolume (particle *particles, double pressure, double boxMatrix[2][2], double *volume, double *totalInteractionEnergy) {
    int result=1;
    double newBoxMatrix[2][2], newTotalInteractionEnergy=0;
    if ((*volume)/VcpPerParticle/N<pressureRealOfNotFluid) {
    //if (pressure>pressureRealOfNotFluid) {  //dozwolone zmiany ksztaltu pudla (faza stala)
        newBoxMatrix[0][0]=boxMatrix[0][0]+(MTRandom0to1(randomStartStep)-0.5)*deltaV;
        newBoxMatrix[1][1]=boxMatrix[1][1]+(MTRandom0to1(randomStartStep)-0.5)*deltaV;
        newBoxMatrix[0][1]=boxMatrix[0][1]+(MTRandom0to1(randomStartStep)-0.5)*deltaV;  //Rectangular volume moves (not romboidal) /**/ around '+(MTRandom0to1(randomStartStep)-0.5)*deltaV'
    } else {    //NIEdozwolone zmiany ksztaltu pudla (faza plynu)
        newBoxMatrix[0][0]=boxMatrix[0][0]+(MTRandom0to1(randomStartStep)-0.5)*deltaV;
        double modifier=newBoxMatrix[0][0]/boxMatrix[0][0];
        newBoxMatrix[1][1]=boxMatrix[1][1]*modifier;
        newBoxMatrix[0][1]=boxMatrix[0][1]*modifier;
    }
    newBoxMatrix[1][0]=newBoxMatrix[0][1];
    double newDetBoxMatrix=newBoxMatrix[0][0]*newBoxMatrix[1][1]-newBoxMatrix[1][0]*newBoxMatrix[0][1],
           newVolume=fabs(newDetBoxMatrix);

    particle particlesInNewBox[activeN];
    cloneParticlesForSpecificBoxMatrix(particlesInNewBox,particles,newBoxMatrix);
    for (int i=0;i<activeN;i++) for (int j=0;j<particles[i].neighCounter;j++) {
        if (i<particles[i].neighbours[j]) {
            getParticlesDistanceSquared(&particlesInNewBox[i],&particlesInNewBox[particles[i].neighbours[j]],newBoxMatrix);
            double drRoot=sqrt(dr[2]);
            if (drRoot<maxDistanceCore) {
                double gamma=atan(dr[1]/dr[0]),
                       aAngle=particles[particles[i].neighbours[j]].phi-gamma,
                       bAngle=particles[i].phi-gamma;
                if (multimerN%2!=0) {
                    //rozważanie, która molekuła jest 'po lewej', a która 'po prawej'
                    if (dr[0]>0) bAngle-=C;
                    else aAngle-=C;
                }
                newTotalInteractionEnergy=checkMoleculeInteriorsOverlap(drRoot,aAngle,bAngle);
                if (newTotalInteractionEnergy<0) {
                    result=0;
                    i=activeN; break;
                }
            }
        }
    }
    if (result) {
        for (int i=0;i<activeN;i++) for (int j=0;j<particles[i].neighCounter;j++) {
            if (i<particles[i].neighbours[j]) {
                getParticlesDistanceSquared(&particlesInNewBox[i],&particlesInNewBox[particles[i].neighbours[j]],newBoxMatrix);
                particlesInNewBox[i].neighEnergyFactor[j]=getEnergy(&particlesInNewBox[i],&particlesInNewBox[particles[i].neighbours[j]]);
                newTotalInteractionEnergy+=particlesInNewBox[i].neighEnergyFactor[j];
            }
        }
        double arg=-((newTotalInteractionEnergy-(*totalInteractionEnergy))/TReduced+pressure*(newVolume-(*volume))-(((double)N+1.0)*log(newVolume/(*volume))+log((newBoxMatrix[0][0]+newBoxMatrix[1][1])/(boxMatrix[0][0]+boxMatrix[1][1]))));  //czynnik ciśnieniowy został podzielony przez kT[=T**eps] (czy raczej: nie został przezeń pomnożony) przy deklaracji ciśnienia, aby tutaj nie dzielić przez kT[=T**eps] przy każdej iteracji
        if (MTRandom0to1(randomStartStep)>exp(arg)) result=0;
        else {
            *volume=newVolume;
            *totalInteractionEnergy=newTotalInteractionEnergy;
            for (int i=0;i<2;i++) for (int j=0;j<2;j++) boxMatrix[i][j]=newBoxMatrix[i][j];
            for (int i=0;i<activeN;i++) {
                for (int j=0;j<2;j++) {
                    particles[i].r[j]=particlesInNewBox[i].r[j];
                    for (int k=0;k<multimerN;k++) particles[i].atomsR[k][j]=particlesInNewBox[i].atomsR[k][j];
                }
                for (int j=0;j<particles[i].neighCounter;j++) if (i<particles[i].neighbours[j]) {
                    particles[i].neighEnergyFactor[j]=particlesInNewBox[i].neighEnergyFactor[j];
                    particles[particles[i].neighbours[j]].neighEnergyFactor[getIndexFromList(particles[particles[i].neighbours[j]].neighbours,particles[particles[i].neighbours[j]].neighCounter,i)]=particlesInNewBox[i].neighEnergyFactor[j];
                }
            }
            detBoxMatrix=newDetBoxMatrix;
        }
    }
    return result;
}
/////////////////  } PARTICLE functions

int createIterationTable () {
    char startArguments[50]; FILE *fileStartArguments = fopen("startArguments.txt","rt");
    if (fileStartArguments==NULL) {printf("Missing file: startArguments.txt\n"); return 1;}
    while (fgets(startArguments,50,fileStartArguments)!=NULL) {
        sscanf(startArguments,"%c",startArguments); char *pEnd;
        iterationTable[fileIterateIterationsNumber][0]=strtod(startArguments,&pEnd);
        iterationTable[fileIterateIterationsNumber][1]=strtod(pEnd,&pEnd);
        iterationTable[fileIterateIterationsNumber++][2]=strtod(pEnd,NULL);
    }
    fclose(fileStartArguments);
    return 0;
}

void addAppendix (char *fileName, char *JOBID, bool jobIdOn) {
    strcpy(buffer,"2D_N-"); strncat(buffer,bufferN,20);
    strncat(buffer,"_gaps-",10); strncat(buffer,bufferGaps,20);
    strncat(buffer,"_G-",5); strncat(buffer,bufferG,5);
    strncat(buffer,"_badanie-",10); strncat(buffer,bufferFolderIndex,5);
    strncat(buffer,"_mN-",5); strncat(buffer,bufferMN,20);
    strncat(buffer,"_mS-",5); strncat(buffer,bufferMS,20);
    strncat(buffer,"_mD-",5); strncat(buffer,bufferMD,20);
    strncat(buffer,"_T-",5); strncat(buffer,bufferTReduced,20);
    mkdir(buffer,S_IRWXU);
    strncat(buffer,"/",2);
    if (jobIdOn) {
        strncat(buffer,JOBID,50);
        strncat(buffer,"_",2);
    }
    strncat(buffer,fileName,200);
    strcpy(fileName,buffer);
}

void adjustOrientationsFile (FILE *file, char *path) {
    if (file==NULL) {
        file=fopen(path,"a"); fprintf(file,"{"); fclose(file);
    } else if (!onlyMath[0]) {
        char bufferForEraseLastChar[200],linia[110]; strcpy(bufferForEraseLastChar,path); strncat(bufferForEraseLastChar,"_BUFF",6);
        FILE *bFELC = fopen(bufferForEraseLastChar,"w");
        int poziomNawiasu=0;
        while (fgets(linia,100,file)!=NULL) {
            sscanf(linia,"%c",linia); int lastIndex=100;
            for (int i=0;i<100;i++) {
                if (linia[i]=='{') poziomNawiasu++;
                else if (linia[i]=='}') poziomNawiasu--;
                if (poziomNawiasu==0) {
                    lastIndex=i;
                    break;
                }
            }
            if (lastIndex==100) fprintf(bFELC,"%s",linia);
            else {
                for (int i=0;i<lastIndex;i++) fprintf(bFELC,"%c",linia[i]);
                fprintf(bFELC,"%c",',');
            }
        }
        fclose(bFELC); fclose(file);
        remove(path); rename(bufferForEraseLastChar,path);
    }
}

void getNextArgument (double prevArg[2], bool countIterations) {
    if (countIterations) if (--iterationsNumber==0) growing=-1;
    if (useFileToIterate) {
        if (++actIteration<fileIterateIterationsNumber) {
            for (int i=0;i<2;i++) prevArg[i]=growing?iterationTable[actIteration][i+1]:iterationTable[fileIterateIterationsNumber-1-actIteration][i+1];
            if (growing) startMinPacFrac=iterationTable[actIteration][0];
            else startMaxPacFrac=iterationTable[fileIterateIterationsNumber-1-actIteration][0];
        } else growing=-1;
    } else if (growing==1) {
        if (multiplyArgument) prevArg[0]*=multiplyFactor;
        else for (int i=0;i<10;i++)
            if (round(prevArg[0]*10000)>=round(intervalMin[i]*10000) && round(prevArg[0]*10000)<round(intervalMax[i]*10000)) {
                double newArg=round(prevArg[0]*10000)+round(intervalDelta[i]*10000);
                prevArg[0]=newArg/10000.0;
                break;
            }
        if (round(prevArg[0]*10000)>round(maxArg*10000)) growing=-1;
    } else if (growing==0) {
        if (multiplyArgument) prevArg[0]/=multiplyFactor;
        else for (int i=0;i<10;i++)
            if (round(prevArg[0]*10000)>round(intervalMin[i]*10000) && round(prevArg[0]*10000)<=round(intervalMax[i]*10000)) {
                double newArg=round(prevArg[0]*10000)-round(intervalDelta[i]*10000);
                prevArg[0]=newArg/10000.0;
                break;
            }
        if (round(prevArg[0]*10000)<round(minArg*10000)) growing=-1;
    }
}

bool isLineCorrect(char linia[4096]) {
    sscanf(linia,"%c",linia);
    int actIndex=0, dataIndex=0; while (dataIndex<3) {
        char data[50]="";
        int licznik=0, dotCounter=0;
        while (linia[actIndex]!='\t' && (dotCounter=linia[actIndex]=='.'?dotCounter+1:dotCounter)<=1 && licznik<50) data[licznik++]=linia[actIndex++];
        if (dotCounter>1 || licznik>=50) {dataIndex=10; continue;} //test of single dot in a number
        actIndex++; dataIndex++;
        if (dataIndex<10 && ((data[0]!='-' && data[1]!='.') || (data[0]=='-' && data[2]!='.'))) dataIndex=10; //test of dot position after first digit
    } if (dataIndex<10 && linia[actIndex]!='\n') dataIndex=10; //test of max 3 numbers in a row

    return (dataIndex<10);
}

double getAvErrorFromSumEps (double sum, double denominator) {
    return sqrt(sum/denominator);
}

void updateTableAndGetActualMean (double table[100], double & mean, int const & changeIndex, double const & changeValue) {
    mean-=table[changeIndex]*0.01; table[changeIndex]=changeValue; mean+=changeValue*0.01;
}

int main(int argumentsNumber, char **arguments) {
/////////////////////////////////////////////// DANE WEJSCIOWE
    int testValue; do {
        char config[500];
        FILE *fileConfig = fopen("config.txt","rt");
        if (fileConfig==NULL) {
            printf("Missing file: config.txt\n");
            return 0;
        }
        int dataIndex=0,intervalLicznik=0;
        while(fgets(config,500,fileConfig)!=NULL) {
            sscanf(config,"%c",config);
            int actIndex=0,licznik=0;
            char data[20]="";
            while (config[actIndex]!='=') actIndex++;
            actIndex++;
            while (config[actIndex]!=';') data[licznik++]=config[actIndex++]; data[licznik]=' ';
            switch (dataIndex) {
                case 0:testValue=strtol(data,NULL,10);break;
                case 1:N=strtol(data,NULL,10);break;
                case 2:gaps=strtol(data,NULL,10);break;
                case 3:multimerN=strtol(data,NULL,10);break;
                case 4:initMode=strtol(data,NULL,10);break;
                case 5:multimerS=strtod(data,NULL);break;
                case 6:multimerD=strtod(data,NULL);break;
                case 7:TReduced=strtod(data,NULL);break;
                case 8:invPotPower=strtod(data,NULL);break;
                case 9:pressureRealOfNotFluid=strtod(data,NULL);break;
                case 10:growing=strtol(data,NULL,10);break;
                case 11:loadedConfiguration=strtol(data,NULL,10);break;
                case 12:loadedArg=strtod(data,NULL);break;
                case 13:{strcpy(loadedJOBID,"j-"); strncat(loadedJOBID,data,licznik);}break;
                case 14:loadedSetStartGenerator=strtol(data,NULL,10);break;
                case 15:loadedSetGenerator=strtol(data,NULL,10);break;
                case 16:iterationsNumber=strtol(data,NULL,10);break;
                case 17:countCollidingPairs=strtol(data,NULL,10);break;
                case 18:intervalSampling=strtol(data,NULL,10);break;
                case 19:intervalOutput=strtol(data,NULL,10);break;
                case 20:saveConfigurations=strtol(data,NULL,10);break;
                case 21:savedConfigurationsInt=strtol(data,NULL,10);break;
                case 22:ODFLength=strtol(data,NULL,10);break;
                case 23:OCFMode=strtol(data,NULL,10);break;
                case 24:neighUpdatingFrequency=strtol(data,NULL,10);break;
                case 25:intervalOrientations=strtol(data,NULL,10);break;
                case 26:skipFirstIteration=strtol(data,NULL,10);break;
                case 27:useSpecificDirectory=strtol(data,NULL,10);break;
                case 28:autoEqLength=strtol(data,NULL,10);break;
                case 29:cyclesOfEquilibration=strtol(data,NULL,10);break;
                case 30:cyclesOfMeasurement=strtol(data,NULL,10);break;
                case 31:intervalResults=strtol(data,NULL,10);break;
                case 32:maxDeltaR=strtod(data,NULL);break;
                case 33:desiredAcceptanceRatioR=strtod(data,NULL);break;
                case 34:desiredAcceptanceRatioV=strtod(data,NULL);break;
                case 35:useFileToIterate=strtol(data,NULL,10);break;
                case 36:startMinPacFrac=strtod(data,NULL);break;
                case 37:startMaxPacFrac=strtod(data,NULL);break;
                case 38:minArg=strtod(data,NULL);break;
                case 39:maxArg=strtod(data,NULL);break;
                case 40:multiplyArgument=strtol(data,NULL,10);break;
                case 41:multiplyFactor=strtod(data,NULL);break;
                default:
                    switch ((dataIndex-42)%3) {
                        case 0: intervalMin[intervalLicznik++/3]=strtod(data,NULL);break;
                        case 1: intervalMax[intervalLicznik++/3]=strtod(data,NULL);break;
                        case 2: intervalDelta[intervalLicznik++/3]=strtod(data,NULL);break;
                    }
                    break;
            }
            dataIndex++;
        }
        fclose(fileConfig);
    } while (testValue!=12345);

    //zlecanie parametrow z poziomu wiersza polecen:
    char JOBID[50]="j-"; int pointNumber=0;
    if (argumentsNumber==1) {
        strncat(JOBID,"none",50);
        if (useFileToIterate) if(createIterationTable()) return 0;
    } else {
        int correctNumberOfArguments=1;
        switch (strtol(arguments[1],NULL,10)) {
            case 0: //ustaw JOBID
                if (argumentsNumber==3) {
                    strncat(JOBID,arguments[2],50); strcpy(loadedJOBID,JOBID);
                    if (useFileToIterate) if(createIterationTable()) return 0;
                } else correctNumberOfArguments=0; break;
            case 1: //ustaw JOBID, singleRun dla parametrow zadanych bezposrednio
                if (argumentsNumber==14) {
                    useFileToIterate=0;
                    strncat(JOBID,arguments[2],50); strcpy(loadedJOBID,JOBID);
                    if (growing) {
                        startMinPacFrac=strtod(arguments[3],NULL); minArg=strtod(arguments[4],NULL);
                    } else {
                        startMaxPacFrac=strtod(arguments[3],NULL); maxArg=strtod(arguments[4],NULL);
                    }
                    N=strtol(arguments[5],NULL,10);
                    gaps=strtol(arguments[6],NULL,10);
                    multimerS=strtod(arguments[7],NULL);
                    multimerD=strtod(arguments[8],NULL);
                    growing=strtol(arguments[9],NULL,10);
                    iterationsNumber=strtol(arguments[10],NULL,10);
                    useSpecificDirectory=strtol(arguments[11],NULL,10);
                    multimerN=strtol(arguments[12],NULL,10);
                    TReduced=strtod(arguments[13],NULL);
                } else correctNumberOfArguments=0; break;
            case 2: //ustaw JOBID, run z najistotniejszymi parametrami z 'config.txt' nadpisanymi z poziomu wywolania
                if (argumentsNumber==15) {
                    if (useFileToIterate) if(createIterationTable()) return 0;
                    strncat(JOBID,arguments[2],50);
                    N=strtol(arguments[3],NULL,10);
                    gaps=strtol(arguments[4],NULL,10);
                    multimerS=strtod(arguments[5],NULL);
                    multimerD=strtod(arguments[6],NULL);
                    growing=strtol(arguments[7],NULL,10);
                    iterationsNumber=strtol(arguments[8],NULL,10);
                    useSpecificDirectory=strtol(arguments[9],NULL,10);
                    skipFirstIteration=strtol(arguments[10],NULL,10);
                    pointNumber=strtol(arguments[11],NULL,10);
                    generatorStartPoint=strtol(arguments[12],NULL,10);
                    multimerN=strtol(arguments[13],NULL,10);
                    TReduced=strtod(arguments[14],NULL);
                } else correctNumberOfArguments=0; break;
            case 3: //ustaw JOBID, tryb loadowany #1 od zadanego argumentu w odpowiednim folderze i trybie
                if (argumentsNumber==15) {
                    if (useFileToIterate) if(createIterationTable()) return 0;
                    strncat(JOBID,arguments[2],50);
                    strcpy(loadedJOBID,"j-"); strncat(loadedJOBID,arguments[3],50);
                    N=strtol(arguments[4],NULL,10);
                    gaps=strtol(arguments[5],NULL,10);
                    multimerS=strtod(arguments[6],NULL);
                    multimerD=strtod(arguments[7],NULL);
                    growing=strtol(arguments[8],NULL,10);
                    loadedConfiguration=1;
                    loadedArg=strtod(arguments[9],NULL);
                    iterationsNumber=strtol(arguments[10],NULL,10);
                    useSpecificDirectory=strtol(arguments[11],NULL,10);
                    skipFirstIteration=strtol(arguments[12],NULL,10);
                    multimerN=strtol(arguments[13],NULL,10);
                    TReduced=strtod(arguments[14],NULL);
                } else correctNumberOfArguments=0; break;
            case 4: //ustaw JOBID, tryb loadowany #2 od zadanego numeru punktu (0->startArg) w odpowiednim folderze i trybie
                if (argumentsNumber==15) {
                    if (useFileToIterate) if(createIterationTable()) return 0;
                    strncat(JOBID,arguments[2],50);
                    strcpy(loadedJOBID,"j-"); strncat(loadedJOBID,arguments[3],50);
                    N=strtol(arguments[4],NULL,10);
                    gaps=strtol(arguments[5],NULL,10);
                    multimerS=strtod(arguments[6],NULL);
                    multimerD=strtod(arguments[7],NULL);
                    growing=strtol(arguments[8],NULL,10);
                    loadedConfiguration=1; loadType=1;
                    pointNumber=strtol(arguments[9],NULL,10);
                    iterationsNumber=strtol(arguments[10],NULL,10);
                    useSpecificDirectory=strtol(arguments[11],NULL,10);
                    skipFirstIteration=strtol(arguments[12],NULL,10);
                    multimerN=strtol(arguments[13],NULL,10);
                    TReduced=strtod(arguments[14],NULL);
                } else correctNumberOfArguments=0; break;
            case 5: //ustaw JOBID, zrob tryb ONLYMATH, gdzie argument wskazuje ile poczatkowych linii Results ma byc pominietych
                if (argumentsNumber==14) {
                    if (useFileToIterate) if(createIterationTable()) return 0;
                    strncat(JOBID,arguments[2],50); strcpy(loadedJOBID,JOBID);
                    onlyMath[0]=1;
                    onlyMath[1]=strtol(arguments[3],NULL,10);
                    N=strtol(arguments[4],NULL,10);
                    gaps=strtol(arguments[5],NULL,10);
                    multimerS=strtod(arguments[6],NULL);
                    multimerD=strtod(arguments[7],NULL);
                    growing=strtol(arguments[8],NULL,10);
                    loadedConfiguration=0;
                    pointNumber=strtol(arguments[9],NULL,10);
                    iterationsNumber=strtol(arguments[10],NULL,10);
                    useSpecificDirectory=strtol(arguments[11],NULL,10);
                    multimerN=strtol(arguments[12],NULL,10);
                    skipFirstIteration=0;
                    TReduced=strtod(arguments[13],NULL);
                } else correctNumberOfArguments=0; break;
            default: {
                printf("Wrong type of run! (0-6)\n");
                return 0;
            } break;
        }
        if (!correctNumberOfArguments) {
            printf("Wrong number of arguments for this type of run!\n");
            printf("If type of run is '0', next arguments: $JOBID\n");
            printf("If type of run is '1', next arguments: $JOBID, startMinPacFrac, minArg, N, gaps, multimerS, multimerD, growing, iterationsNumber, useSpecificDirectory, multimerN, TReduced\n");
            printf("If type of run is '2', next arguments: $JOBID, N, gaps, multimerS, multimerD, growing, iterationsNumber, useSpecificDirectory, skipFirstIteration, pointNumber, generatorStartPoint, multimerN, TReduced\n");
            printf("If type of run is '3', next arguments: $JOBID, JOBID of configuration to load, N, gaps, multimerS, multimerD, growing, loadedArg, iterationsNumber, useSpecificDirectory, skipFirstIteration, multimerN, TReduced\n");
            printf("If type of run is '4', next arguments: $JOBID, JOBID of configuration to load, N, gaps, multimerS, multimerD, growing, pointNumber, iterationsNumber, useSpecificDirectory, skipFirstIteration, multimerN, TReduced\n");
            printf("If type of run is '5', next arguments: $JOBID, lines to skip from Results, N, gaps, multimerS, multimerD, growing, pointNumber, iterationsNumber, useSpecificDirectory, multimerN, TReduced\n");
            return 0;
        }
    }

    //ostatnie konfiguracje
    if (useFileToIterate) {
        if (growing) {
            startMinPacFrac=iterationTable[0][0];
            for (int i=0;i<2;i++) startArg[i]=iterationTable[0][i+1];
        } else {
            startMaxPacFrac=iterationTable[fileIterateIterationsNumber-1][0];
            for (int i=0;i<2;i++) startArg[i]=iterationTable[fileIterateIterationsNumber-1][i+1];
        }
    } else {
        startArg[0]=growing?minArg:maxArg;
        startArg[1]=0;
    }

    //pressureRealOfNotFluid - dostosowanie do anisotropii, tutaj używane jako graniczne v* (w przypadku używania jako p* należy pamiętać o komendzie <<pressureRealOfNotFluid/=(multimerD*multimerD)>>)
    /*switch (multimerN) {  //zadane będzie pRONF=0.5, dla badań histerez
        case 6: switch ((int)round(multimerD/multimerS*100)) {
            case 50: pressureRealOfNotFluid=1.72; break; //0.5
            case 60: pressureRealOfNotFluid=1.48; break; //0.6
            case 75: pressureRealOfNotFluid=1.41; break; //0.75
            case 100: pressureRealOfNotFluid=1.35; break; //1
            case 150: pressureRealOfNotFluid=1.29; break; //1.5
            case 300: pressureRealOfNotFluid=1.24; break; //3
            case 100000: pressureRealOfNotFluid=1.17; break; //1000
        } break;
    }*/

    deltaR=maxDeltaR*multimerD; deltaV*=multimerD;  //jednostka odległości to d^**\sigma
    for (int i=0;i<pointNumber;i++) getNextArgument(startArg,false);
    if (loadedConfiguration && loadType) loadedArg=startArg[0];
    activeN=N-gaps;

    if (N%56!=0 && N%780!=0 && fabs(sqrt(N)-floor(sqrt(N)))>0.000001) {
        printf("ERROR: Not supported N: %d.\n",N);
        return 0;
    }

    //stale wynikajace z zadanych parametrow multimerow
    L=multimerS/multimerD;
    C=pi/(double)multimerN; deltaPhi=deltaR*2.0*sin(C)/multimerS;
    ROkreguOpisanego=multimerS/(2.0*sin(C));
    hTrojkataWielokata=ROkreguOpisanego*cos(C);
    absoluteMinimum=atan(L/(2.0*L/tan(C)+sqrt(4.0-L*L)));
    absoluteMinimum2=C-absoluteMinimum;
    switch ((int)(multimerD*1000000)) {
        case 500001: minDistance=1.8037364284584916; break;
        case 501000: minDistance=1.8331940077554698; break;
        case 502000: minDistance=1.845827233464671; break;
        case 503000: minDistance=1.8555403664338204; break;
        case 504000: minDistance=1.8637442860393154; break;
        case 505000: minDistance=1.8709851905229433; break;
        case 506000: minDistance=1.8775430588681767; break;
        case 507000: minDistance=1.8835841266820283; break;
        case 508000: minDistance=1.8892166507639028; break;
        case 509000: minDistance=1.8945157885020365; break;
        case 510000: minDistance=1.899536233850406; break;
        case 511000: minDistance=1.9043192554184105; break;
        case 512000: minDistance=1.9088968985500718; break;
        case 513000: minDistance=1.913294634275949; break;
        case 514000: minDistance=1.9175331038931438; break;
        case 515000: minDistance=1.9216293101072457; break;
        case 516000: minDistance=1.9255974552691306; break;
        case 517000: minDistance=1.9294495466426982; break;
        case 518000: minDistance=1.9331958432559069; break;
        case 519000: minDistance=1.9368451922349774; break;
        case 520000: minDistance=1.9404052862930774; break;
        default: if (multimerN%2==0) minDistance=getMinimalDistanceAnalyticalMethodForEvenHCM(absoluteMinimum,absoluteMinimum);
                 else minDistance=getMinimalDistanceAnalyticalMethodForOddHCM(absoluteMinimum,absoluteMinimum-C); break;
    }
    maxDistance=ROkreguOpisanego*2+multimerD; maxDistanceCore=ROkreguOpisanego*2;

    //nazwy folderow na podstawie parametrow programu
    sprintf(bufferG,"%d",growing); sprintf(bufferN,"%d",N); sprintf(bufferGaps,"%d",gaps);
    sprintf(bufferMN,"%d",multimerN); sprintf(bufferMS,"%.2f",multimerS); sprintf(bufferMD,"%.6f",multimerD); sprintf(bufferTReduced,"%.3E",TReduced);

    int folderIndex=useSpecificDirectory, checkNext;
    char bufferCheckFolderExisting[200];
    FILE *checkFolderExisting;
    if (!folderIndex) do {
        sprintf(bufferFolderIndex,"%d",++folderIndex);
        strcpy(bufferCheckFolderExisting,resultsFileName);
        addAppendix(bufferCheckFolderExisting,JOBID,false);
        checkFolderExisting = fopen(bufferCheckFolderExisting,"rt");
        if (checkFolderExisting!=NULL) {
            fclose(checkFolderExisting);
            checkNext=1;
        } else checkNext=0;
    } while (checkNext);
    sprintf(bufferFolderIndex,"%d",folderIndex);
    addAppendix(resultsFileName,JOBID,false);
    addAppendix(excelResultsFileName,JOBID,false);
    addAppendix(configurationsFileName,JOBID,true);
    addAppendix(loadConfigurationsFileName,loadedJOBID,true); strncat(loadConfigurationsFileName,"_arg-",6); sprintf(buffer,"%.4E",loadedArg); strncat(loadConfigurationsFileName,buffer,100); strncat(loadConfigurationsFileName,".txt",5);
    addAppendix(orientationsFileName,JOBID,true);
    if (OCFMode) addAppendix(orientatCorrelFunFileName,JOBID,true);
    addAppendix(orientationsResultsFileName,JOBID,true);
    addAppendix(configurationsListFileName,JOBID,false);

    particle particles[N];
    double args[10];

    FILE *fileResults, *fileExcelResults, *fileConfigurations, *fileSavedConfigurations, *fileOrientations, *fileOrientatCorrelFun, *fileConfigurationsList, *fileAllResults, *fileAllOrientations, *fileOrientationsResults, *fileAllOrientationsResults;
    fileResults = fopen(resultsFileName,"rt"); if (fileResults==NULL) {
        fileResults = fopen(resultsFileName,"a");
        if (saveConfigurations) fprintf(fileResults,"Cycles\tPressure*\tVolume\tdVolume\tBoxMatrix[0][0]\tdBoxMatrix[0][0]\tBoxMatrix[1][1]\tdBoxMatrix[1][1]\tBoxMatrix[1][0]([0][1])\tdBoxMatrix[1][0]([0][1])\tRho\tdRho\tV/V_cp\tdV/V_cp\tS1111\tdS1111\tS1122\tdS1122\tS1212\tdS1212\tS2222\tdS2222\tS1112\tdS1112\tS1222\tdS1222\tavNu\tdAvNu\tavNu2\tdAvNu2\tavB*\tdAvB*\tavMy*\tdAvMy*\tavE*\tdAvE*\tODFMax_One\t<cos(6Phi)>_One\tODFMax_All\t<cos(6Phi)>_All\tPhiOfODFMax_All\tavPhi_All\tdPhiCyclesInterval\tavAbsDPhi\n");
        else fprintf(fileResults,"Cycles\tPressure*\tVolume\tdVolume\tBoxMatrix[0][0]\tdBoxMatrix[0][0]\tBoxMatrix[1][1]\tdBoxMatrix[1][1]\tBoxMatrix[1][0]([0][1])\tdBoxMatrix[1][0]([0][1])\tRho\tdRho\tV/V_cp\tdV/V_cp\tS1111\tdS1111\tS1122\tdS1122\tS1212\tdS1212\tS2222\tdS2222\tS1112\tdS1112\tS1222\tdS1222\tavNu\tdAvNu\tavNu2\tdAvNu2\tavB*\tdAvB*\tavMy*\tdAvMy*\tavE*\tdAvE*\tODFMax_One\t<cos(6Phi)>_One\tODFMax_All\t<cos(6Phi)>_All\tPhiOfODFMax_All\tavPhi_All\n");
        fclose(fileResults);
    }
    fileExcelResults = fopen(excelResultsFileName,"rt"); if (fileExcelResults==NULL) {
        fileExcelResults = fopen(excelResultsFileName,"a");
        if (saveConfigurations) fprintf(fileExcelResults,"Pressure*\tV/V_cp\tavNu\tavNu2\tODFMax_All\t<cos(6Phi)>_All\tPhiOfODFMax_All\tavPhi_All\tavB*\tavMy*\tavE*\tdPhiCyclesInterval\tavAbsDPhi\n");
        else fprintf(fileExcelResults,"Pressure*\tV/V_cp\tavNu\tavNu2\tODFMax_All\t<cos(6Phi)>_All\tPhiOfODFMax_All\tavPhi_All\tavB*\tavMy*\tavE*\n");
        fclose(fileExcelResults);
    }




/////////////////////////////////////////////// WARUNKI POCZATKOWE

    double arg[2]={startArg[0],startArg[1]}, oldBoxMatrix[2][2], oldTotalInteractionEnergy;
    while (growing>=0) {
        double totalInteractionEnergy, pressureReduced=arg[0], pressure=pressureReduced/multimerD/multimerS, //pressureReduced=\tau*d^**\sigma^2/kT - Dwa punkty widzenia: 1) zmniejszenie sigma ZWIĘKSZA JEDNOSTKE p^*, zatem ten sam STAN FIZYCZNY jak przy sigma=1 bedzie przy mniejszym p^*. pReal to tak naprawde pReducedAdjusted. Inny, równoważny punkt widzenia, to 2) pReal redukuje objetosc, ktora NIE jest wyrazana w jednostkach sigma. Objetosc jest obliczana z boxMatrix, ktory jest inicjowany z czynnikiem *sigma, zatem MA jednostkę (bo niestety sigma jest zdefiniowana jako zmienna i DA się ją zmieniać), a NIE jest zredukowany. Przeciez gdyby sigma=2, to boxMatrix bylby 2x wiekszy, a 'w jednostkach sigma' (zredukowany) powinien pozostac identyczny
               boxMatrix[2][2],matrixOfParticlesSize[2],unitCellAtCP[2],                                //obydwa sprowadzają się do tego, że przy liczeniu prawdopodobieństwa ma być jednostka zredukowana (bezwymiarowa): 1) zakłada, że volume jest zredukowane, więc dostosowuje pReduced do stanu fizycznego; 2) zakłada, że pReduced już jest OK (w końcu jest zredukowane), tylko po prostu objętość NIE jest zredukowana, i trzeba ją zredukować dzieląc przez sigma^2
               matrixCellXY[6][6][2],matrixCellPhi[6][6];                                               //kwestia kT[=T**eps]: \tau (czyli zmienna pressure) przy deklaracji jest dzielona przez kT[=T**eps] (czy raczej: po prostu nie jest przezeń mnożona), aby przy obliczeniach zmiany entalpii nie dzielić przez kT[=T**eps] "co iterację"
        int n[2]; //n[X/Y], matrixCell[n[XMax]][n[YMax]][x/y], zatem: n[X/Y](max)=6                     //przy układach twardych czynnik kT nie ma znaczenia [jeżeli pressure zdefiniuje się uwzględniając *kT, to przy obliczeniach entalpii się przez niego dzieli], natomiast w miękkich ma, bo czynnik kT wchodzi w mianownik energii potencjalnej [ktora przy twardych układach jest 0 lub \infty]. Stąd, w twardych układach naturalne jest branie p^*=p*sig^2/kT, natomiast w miękkich już raczej coś w stylu: p^*=p*sig^2/eps [definiuje się jednostkę energii =eps].
        switch (multimerN) {                                                                            //To jednak tylko numeryka, można (dla wygody) dalej definiować p^* w oparciu o kT - wówczas zmiana T* symulacji sprawi, że przy jednakowym zestawie p^* (startArguments.txt) badane będą różne p (co może być porządane albo nie). Bez zdefiniowania jednostki energii nie da się jednak sensownie wyrazić T^* (chyba, że np. jednostkę temperatury się wprowadzi jakąś t0 i T^*=T/t0, ale tak się raczej nie robi). Podstawowymi jednostkami, na których wyraża
            case 6: {//dla heksamerow o dowolnym d/\sigma                                               //się zredukowane są: długość[sigma], energia[eps] i masa[w tym programie nie ma potrzeby jej definiować]. UWAGA: jednostki oczywiście działają 'w domyśle' (jako 1). NIE powinno być ich w programie w postaci zmiennych (istnienie zmiennej sigma czy multimerS to tutaj relikt błędu przeszłości, 'eps' już nie ma swojej zmiennej). Jeżeli byłaby później potrzeba: zawsze można wyniki (EOS) przetransformować do postaci p^*=p*sig^2/eps (=p*sig^2*T^*/kT)
                unitCellAtCP[0]=minDistance; unitCellAtCP[1]=sqrt(3)*minDistance;
                n[0]=1; n[1]=2;
                matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=absoluteMinimum2;
                matrixCellXY[0][1][0]=unitCellAtCP[0]/2.0; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=absoluteMinimum2;
            } break;
            case 5: {//dla pentamerow o d/\sigma=1
                unitCellAtCP[0]=2.4048671732*multimerS; unitCellAtCP[1]=4.2360679772*multimerS;
                n[0]=1; n[1]=2;
                matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=C;
                matrixCellXY[0][1][0]=1.0131106571*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=0;
            } break;
            case 7: switch (initMode) {//dla heptamerow o d/\sigma=1, struktury jak w WojTreKow2003PRE
                case 0: {//struktura A
                    unitCellAtCP[0]=3.0566685376*multimerS; unitCellAtCP[1]=5.5150210832*multimerS;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=C;
                    matrixCellXY[0][1][0]=1.9144263193*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=0;
                } break;
                case 1: {//struktura B
                    unitCellAtCP[0]=3.0566685376*multimerS; unitCellAtCP[1]=5.5382990167*multimerS;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=C;
                    matrixCellXY[0][1][0]=1.3656999121*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=0;
                } break;
                case 2: {//struktura C
                    unitCellAtCP[0]=6.0309063912*multimerS; unitCellAtCP[1]=5.5371391576*multimerS;
                    n[0]=2; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=-0.1644029030;
                    matrixCellXY[1][0][0]=unitCellAtCP[0]/2.0; matrixCellXY[1][0][1]=0.1230590461*multimerS; matrixCellPhi[1][0]=0.1644029030;
                    matrixCellXY[0][1][0]=1.2832327269*multimerS; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=-0.2843960475;
                    matrixCellXY[1][1][0]=4.2986859225*multimerS; matrixCellXY[1][1][1]=2.8916286250*multimerS; matrixCellPhi[1][1]=0.2843960475;
                } break;
            } break;
            case 3: switch (initMode) {//struktury trimerow jak z pracy KVT
                case 0: {//struktura INC (Isotropic Nonchiral Crystal)
                    unitCellAtCP[0]=0.5*multimerS*sqrt(3)+sqrt(multimerD*multimerD-multimerS*multimerS*0.25); unitCellAtCP[1]=unitCellAtCP[0]*sqrt(3);
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=C;
                    matrixCellXY[0][1][0]=unitCellAtCP[0]/2.0; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=C;
                } break;
                case 1: {//struktura R2 (Rectangular Crystal with 2 molecules in the unit cell)
                    unitCellAtCP[0]=0.5*multimerS*sqrt(3)+sqrt(multimerD*multimerD-multimerS*multimerS*0.25); unitCellAtCP[1]=2.0*cos(C/2.0-absoluteMinimum2)*minDistance;
                    n[0]=1; n[1]=2;
                    matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=0;
                    matrixCellXY[0][1][0]=-sin(C/2.0-absoluteMinimum2)*minDistance; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=C;
                } break;
            } break;
            default: {
                unitCellAtCP[0]=maxDistance; unitCellAtCP[1]=maxDistance*sqrt(3);
                n[0]=1; n[1]=2;
                matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0; matrixCellPhi[0][0]=C;
                matrixCellXY[0][1][0]=unitCellAtCP[0]/2.0; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0; matrixCellPhi[0][1]=C;
            } break;
        } VcpPerParticle=unitCellAtCP[0]*unitCellAtCP[1]/(double)n[0]/(double)n[1];
        if (N%56==0) {matrixOfParticlesSize[0]=7; matrixOfParticlesSize[1]=8;}
        else if (N%780==0) {matrixOfParticlesSize[0]=26; matrixOfParticlesSize[1]=30;}
        else if (floor(sqrt(N))==sqrt(N)) {matrixOfParticlesSize[0]=matrixOfParticlesSize[1]=sqrt(N);}
        double NLinearMod = sqrt(N/matrixOfParticlesSize[0]/matrixOfParticlesSize[1]);
        if (startArg[0]==arg[0]) {
            for (int i=0;i<2;i++) boxMatrix[i][i]=matrixOfParticlesSize[i]*unitCellAtCP[i]/(double)n[i]*NLinearMod*(growing?sqrt(startMinPacFrac):sqrt(startMaxPacFrac));
            boxMatrix[1][0]=0.0; boxMatrix[0][1]=0.0;
        } else {
            for (int i=0;i<2;i++) for (int j=0;j<2;j++) boxMatrix[i][j]=oldBoxMatrix[i][j];
            totalInteractionEnergy=oldTotalInteractionEnergy;
        }
        detBoxMatrix=boxMatrix[0][0]*boxMatrix[1][1]-boxMatrix[1][0]*boxMatrix[0][1];
        double volume=fabs(detBoxMatrix), rho=N/volume, pacFrac=1.0/VcpPerParticle/rho;

        if (!onlyMath[0]) {
            if (arg[0]==startArg[0] && !loadedConfiguration) {
                printf("INIT POS.- N: %d, gaps: %d, growing: %d, StartPressRed: %.7E (StartDen: %.7E, startPacFrac: %.7E, T*: %.3E, invPotPower: %.2f), mN: %d, mS: %.2f, mD: %.6f\n",N,gaps,growing,startArg[0],rho,pacFrac,TReduced,invPotPower,multimerN,multimerS,multimerD);
                if (!initPositions(particles,boxMatrix,matrixOfParticlesSize,n,matrixCellXY,matrixCellPhi,pacFrac,volume)) return 0;
                adjustAngles(particles,boxMatrix);   //dla układów, w których obracanie cząstek pomaga uzyskać initowy układ (np. HCH)
                //adjustOddRowsTranslation(particles,boxMatrix,(int)(matrixOfParticlesSize[0]*NLinearMod));   //dla układów, w których translacja co drugiego rzędu pomaga uzyskać initowy układ (np. HCT)
                updateNeighbourList(particles,boxMatrix,volume);
                totalInteractionEnergy=getEnergyAll(particles,boxMatrix,2);
            } else if (loadedConfiguration) {
                char configurations[4096];
                FILE *fileCTL = fopen(loadConfigurationsFileName,"rt");
                if (fileCTL==NULL) {
                    printf("Missing file (configuration): %s\n",loadConfigurationsFileName);
                    return 0;
                }
                for (int i=0;i<3;i++) fgets(configurations,4096,fileCTL);
                int character,dataType=0,pIndex=-1; char data[50]=""; int actIndex=0;
                while (dataType<11) {
                    character=fgetc(fileCTL); //character is in int, but it can work as char
                    if (dataType<10) { //stage #1 configuration parameters
                        if (character==' ') {
                            data[actIndex++]=' '; //end of data without clearing the entire array
                            args[dataType++]=strtod(data,NULL);
                            if (dataType==10) {
                                boxMatrix[0][0]=args[5]; boxMatrix[1][1]=args[6]; boxMatrix[1][0]=args[7]; boxMatrix[0][1]=args[7];
                                deltaR=args[8]; deltaPhi=deltaR*2.0*sin(C)/multimerS; deltaV=args[9];
                                arg[0]=args[3]; pressureReduced=arg[0]; pressure=pressureReduced/multimerD/multimerS;

                                detBoxMatrix=boxMatrix[0][0]*boxMatrix[1][1]-boxMatrix[1][0]*boxMatrix[0][1];
                                volume=fabs(detBoxMatrix); rho=N/volume; pacFrac=1.0/VcpPerParticle/rho;
                            }
                            actIndex=0;
                            continue;
                        }
                        data[actIndex++]=character;
                    } else { //stage #2 configuration parameters (coordinates of particles)
                        if (character=='m') {pIndex++; continue;}
                        if (pIndex>=0) {
                            for (int i=0;i<3;i++) {
                                character=fgetc(fileCTL);
                                while (character!=',') {
                                    data[actIndex++]=character;
                                    character=fgetc(fileCTL);
                                } data[actIndex++]=' '; actIndex=0;
                                if (i<2) particles[pIndex].r[i]=strtod(data,NULL);
                                else particles[pIndex].phi=strtod(data,NULL);
                            }
                            particles[pIndex].normR[0]=(boxMatrix[1][1]*particles[pIndex].r[0]-boxMatrix[0][1]*particles[pIndex].r[1])/detBoxMatrix;
                            particles[pIndex].normR[1]=-(boxMatrix[1][0]*particles[pIndex].r[0]-boxMatrix[0][0]*particles[pIndex].r[1])/detBoxMatrix;
                            computeAtomsPositions(&particles[pIndex]);
                            while (character!=']') character=fgetc(fileCTL); fgetc(fileCTL); //next to read: 'm'
                            if (pIndex>=activeN-1) dataType++;
                        }
                    }
                }
                fclose(fileCTL);
                printf("LOADING POS.- N: %d, gaps: %d, growing: %d, StartPressRed: %.7E (startDen: %.7E, startPacFrac: %.7E, T*: %.3E, invPotPower: %.2f), RandStart: %.1f, RandStep: %.1f, Cycles: %ld, DeltaR: %.4E, DeltaV: %.4E\n",N,gaps,growing,args[3],args[2],pacFrac,TReduced,invPotPower,args[0],args[1],(long)args[4],args[8],args[9]);
                //for (int i=0;i<2;i++) for (int j=0;j<2;j++) printf("boxMatrix[%d][%d]=%.17E\n",i,j,boxMatrix[i][j]);
                //for (int i=0;i<activeN;i++) printf("%d: %.17E,  %.17E,  %.17E\n",i,particles[i].r[0],particles[i].r[1],particles[i].phi);return 0;
                updateNeighbourList(particles,boxMatrix,volume);
                totalInteractionEnergy=getEnergyAll(particles,boxMatrix,0);
                printf("Checking molecule's core overlaps in loaded file... "); fflush(stdout);
                if (totalInteractionEnergy<0) {
                    printf("Configuration from loaded file contains molecule's core overlap(s) [energy->infinity].\n");
                    return 0;
                } else  {printf("done\n"); fflush(stdout);}
            }
        }

        if (skipFirstIteration) {
            printf("Skipping first iteration...\n");
            skipFirstIteration=0;
        } else {
            if (!onlyMath[0]) {
                if (loadedConfiguration) {
                    if (loadedSetStartGenerator) {
                        printf("Setting start position of p-random number generator to position from file...\n");
                        InitMT((unsigned int)args[0]);
                        randomStartStep[0]=args[0];
                    } else {
                        printf("Setting start position of p-random number generator to position from file - DISABLED\n");
                        if (generatorStartPoint==0) {
                            generatorStartPoint=time(0);
                            printf("Setting start position of p-random number generator to actual CPU time...\n");
                        } else printf("Setting start position of p-random number generator to %ld...\n",generatorStartPoint);
                        InitMT((unsigned int)generatorStartPoint);
                        randomStartStep[0]=generatorStartPoint;
                    }
                    randomStartStep[1]=0;
                    if (loadedSetGenerator) {
                        printf("Setting p-random number generator to last position from file...\n");
                        for (double i=0;i<args[1];i++) MTGenerate(randomStartStep);
                    } else printf("Setting p-random number generator to last position from file - DISABLED\n");
                } else {
                    if (generatorStartPoint==0) {
                        generatorStartPoint=time(0);
                        printf("Setting start position of p-random number generator to actual CPU time...\n");
                    } else printf("Setting start position of p-random number generator to %ld...\n",generatorStartPoint);
                    InitMT((unsigned int)generatorStartPoint);
                    randomStartStep[0]=generatorStartPoint;
                    randomStartStep[1]=0;
                }
                printf("Start of equilibration at reduced pressure: %.7E (startDen: %.7E, startPacFrac: %.7E)...",pressureReduced,rho,pacFrac);
                if (autoEqLength==0) printf(" (%ld cycles)\n",cyclesOfEquilibration); else printf(" (auto, min. length: %ld cycles)\n",cyclesOfEquilibration*intervalSampling);
            } else printf("Start of mathOnly mode for: N: %d, gaps: %d, growing: %d, pressRed: %.7E, mN: %d, mS: %.2f, mD: %.6f\n",N,gaps,growing,pressureReduced,multimerN,multimerS,multimerD);
            fflush(stdout);




/////////////////////////////////////////////// RDZEN MC

            long volumeMoveChance=(int)ceil(activeN/sqrt(activeN)),
                fullCycle=activeN+volumeMoveChance,cycle=0,   //UWAGA cycle LONG, nie moze byc za duzo cykli
                timeStart,timeEquilibration=0,timeEnd,
                attemptedNumberR=0, displacedNumberR=0,
                attemptedNumberV=0, displacedNumberV=0,
                cyclesOfMeasurementBuffer=arg[1]==0?cyclesOfMeasurement:0,
                cyclesOfEquilibrationBuffer=cyclesOfEquilibration;
            double deltaRTable[100], deltaRMean=deltaR, deltaVTable[100], deltaVMean=deltaV;
            for (int i=0;i<100;i++) {deltaRTable[i]=deltaRMean; deltaVTable[i]=deltaVMean;}
            int simulationStage=cyclesOfEquilibrationBuffer>0?0:cyclesOfMeasurementBuffer>0?1:2;  //0-equilibration, 1-measurement, 2-end
            double autoEqVolumeTable[cyclesOfEquilibrationBuffer],autoEqBalance=0;
            autoEqVolumeTable[0]=volume; for (int i=1;i<cyclesOfEquilibrationBuffer;i++) autoEqVolumeTable[i]=-1;
            int volumeMove=0, cycleCounter=0, indexScanned=(matrixOfParticlesSize[0]*round(matrixOfParticlesSize[1]*NLinearMod/2.0)-round(matrixOfParticlesSize[0]/2.0))*NLinearMod,autoEqCounter=0,autoEqCheck=0;

            char allResultsFileName[200],bufferConfigurations[200],bufferSavedConfigurations[200],bufferOrientations[200],bufferOrientatCorrelFun[200],allOrientationsFileName[200],bufferOrientationsResults[200],allOrientationsResultsFileName[200],bufferPressure[100];
            strcpy(allResultsFileName,configurationsFileName); strcpy(bufferConfigurations,configurationsFileName);
            strcpy(bufferOrientations,orientationsFileName); strcpy(allOrientationsFileName,orientationsFileName);
            if (OCFMode) strcpy(bufferOrientatCorrelFun,orientatCorrelFunFileName);
            strcpy(bufferOrientationsResults,orientationsResultsFileName); strcpy(allOrientationsResultsFileName,orientationsResultsFileName);
            sprintf(bufferPressure,"%.4E",pressureReduced);
            strncat(allResultsFileName,"_arg-",6); strncat(allResultsFileName,bufferPressure,100); strncat(allResultsFileName,"_Results.txt",13);
            strncat(bufferConfigurations,"_arg-",6); strncat(bufferConfigurations,bufferPressure,100); strcpy(bufferSavedConfigurations,bufferConfigurations); strncat(bufferConfigurations,".txt",5); strncat(bufferSavedConfigurations,"_transient.txt",15);
            strncat(bufferOrientations,"_arg-",6); strncat(bufferOrientations,bufferPressure,100); strncat(bufferOrientations,".txt",5);
            strncat(allOrientationsFileName,"_arg-",6); strncat(allOrientationsFileName,bufferPressure,100); strncat(allOrientationsFileName,"_allOnt.txt",12);
            if (OCFMode) {
                strncat(bufferOrientatCorrelFun,"_arg-",6); strncat(bufferOrientatCorrelFun,bufferPressure,100); strncat(bufferOrientatCorrelFun,".txt",5);
            }
            strncat(bufferOrientationsResults,"_arg-",6); strncat(bufferOrientationsResults,bufferPressure,100); strncat(bufferOrientationsResults,".txt",5);
            strncat(allOrientationsResultsFileName,"_arg-",6); strncat(allOrientationsResultsFileName,bufferPressure,100); strncat(allOrientationsResultsFileName,"_allOnt.txt",12);

            fileAllResults = fopen(allResultsFileName,"a");
            adjustOrientationsFile(fileOrientations=fopen(bufferOrientations,"rt"),bufferOrientations);
            fileOrientations = fopen(bufferOrientations,"a");
            fileAllOrientations = fopen(allOrientationsFileName,"a");
            if (saveConfigurations) fileSavedConfigurations = fopen(bufferSavedConfigurations,"a");
            if (OCFMode) {
                adjustOrientationsFile(fileOrientatCorrelFun=fopen(bufferOrientatCorrelFun,"rt"),bufferOrientatCorrelFun);
                fileOrientatCorrelFun = fopen(bufferOrientatCorrelFun,"a");
            }
            if (onlyMath[0]) simulationStage=2;

            timeStart=time(0);
            while (simulationStage<2) {
                int randIndex;
                if (volumeMove) {
                    randIndex = (int)(MTRandom0to1(randomStartStep)*activeN);
                    volumeMove=0;
                } else randIndex = (int)(MTRandom0to1(randomStartStep)*fullCycle);
                if (randIndex<activeN) {
                    attemptedNumberR++;
                    if (attemptToDisplaceAParticle(particles,randIndex,boxMatrix,&totalInteractionEnergy))
                        displacedNumberR++;
                } else {
                    volumeMove=1;
                    attemptedNumberV++;
                    if (attemptToChangeVolume(particles,pressure,boxMatrix,&volume,&totalInteractionEnergy))
                        displacedNumberV++;
                }

                cycleCounter++;
                if (cycleCounter>=fullCycle) {
                    cycleCounter=0;
                    cycle++;

                    if (cycle%intervalSampling==0) {
                        if (simulationStage==0) {
                            if (autoEqLength==0) {
                                if (cycle>cyclesOfEquilibrationBuffer) simulationStage=1;
                            } else {  //metoda 'przyrostow' w sposob 'listy' przy uzyciu tablicy-ale bez przepisywania wszystkich elementow (o 1 w dol) jezeli tablica sie zapelnia
                                autoEqBalance+=volume-autoEqVolumeTable[autoEqCounter++];
                                if (autoEqCounter>=cyclesOfEquilibrationBuffer) autoEqCounter=0;
                                int nextIndex=autoEqCounter+1; if (nextIndex>=cyclesOfEquilibrationBuffer) {nextIndex=0; if (autoEqCheck==0) autoEqCheck=autoEqBalance>0?1:-1;}
                                if (autoEqVolumeTable[autoEqCounter]>0) autoEqBalance-=autoEqVolumeTable[nextIndex]-autoEqVolumeTable[autoEqCounter];
                                autoEqVolumeTable[autoEqCounter]=volume;
                                if (autoEqCheck!=0 && autoEqBalance*autoEqCheck<0) {
                                    simulationStage=1; cyclesOfEquilibrationBuffer=cycle;
                                }
                            }
                            if (simulationStage==1) {
                                printf("Equilibration finished after: %ld cycles (%ldsec).\n",cyclesOfEquilibrationBuffer,time(0)-timeStart);
                                fflush(stdout);
                            }
                        }
                        double acceptanceRatioR = displacedNumberR/(double)attemptedNumberR,
                               acceptanceRatioV = displacedNumberV/(double)attemptedNumberV;
                        if (cycle%neighUpdatingFrequency==0) {
                            updateNeighbourList(particles,boxMatrix,volume);
                            totalInteractionEnergy=getEnergyAll(particles,boxMatrix,2);
                        }

                        /////wypisywanie danych czesciej niz normalnie i PRZED zrownowagowaniem
                        /*if (cycle%50==0) {
                            int collidingPairs=0;                         
                            if (countCollidingPairs) {
                                for (int i=0;i<activeN-1;i++) for (int j=i+1;j<activeN;j++) {
                                    getParticlesDistance(particles[i],particles[j],boxMatrix);
                                    double drRoot=sqrt(dr[2]);
                                    if (drRoot<maxDistance) {
                                        double gamma=atan(dr[1]/dr[0]),
                                               aAngle=particles[j].phi-gamma,
                                               bAngle=particles[i].phi-gamma;
                                        if (multimerN%2!=0) {
                                            //rozważanie, która molekuła jest 'po lewej', a która 'po prawej'
                                            if (dr[0]>0) bAngle-=C;
                                            else aAngle-=C;
                                        }
                                        if (checkMoleculeInteriorsOverlap(drRoot,aAngle,bAngle)<0) collidingPairs++;
                                    }
                                }
                            }

                            rho=N/volume; pacFrac=1.0/VcpPerParticle/rho;
                            if (countCollidingPairs) printf("Cycle: %ld, CollPairs: %d\n",(cycle+args[4]),collidingPairs);
                            else printf("Cycle: %ld\n",(cycle+args[4]));
                            printf("   AccRatR: %.4E, dR: %.4E, AccRatV: %.4E, dV: %.4E\n",acceptanceRatioR,deltaR,acceptanceRatioV,deltaV);
                            printf("   Dens: %.4E, V/V_cp: %.4E, PressRed: %.4E\n",rho,pacFrac,pressureReduced);
                            printf("   box00: %.8E, box11: %.8E, box01(10): %.8E\n",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1]);
                            fflush(stdout);
                        }*/
                        /////

                        if (simulationStage==1) {
                            if (timeEquilibration==0) {
                                timeEquilibration=time(0);

                                printf("Checking molecule's core overlaps in equilibrated configuration... "); fflush(stdout);
                                if (getEnergyAll(particles,boxMatrix,1)<0) {
                                    printf("Equilibrated configuration contains molecule's core overlap(s) [energy->infinity]. ");
                                    char allResultsErrorFileName[200];
                                    strcpy(allResultsErrorFileName,allResultsFileName); strncat(allResultsErrorFileName,".err",5);
                                    if (rename(allResultsFileName,allResultsErrorFileName)==0) printf("Results file successfully renamed (.err).\n");
                                    else printf("Error renaming results file (.err).\n");
                                    return 0;
                                } else  {printf("done\n"); fflush(stdout);}
                            }

                            int collidingPairs=0;
                            if (countCollidingPairs) {
                                for (int i=0;i<activeN-1;i++) for (int j=i+1;j<activeN;j++) {
                                    getParticlesDistanceSquared(&particles[i],&particles[j],boxMatrix);
                                    double drRoot=sqrt(dr[2]);
                                    if (drRoot<maxDistanceCore) {
                                        double gamma=atan(dr[1]/dr[0]),
                                               aAngle=particles[j].phi-gamma,
                                               bAngle=particles[i].phi-gamma;
                                        if (multimerN%2!=0) {
                                            //rozważanie, która molekuła jest 'po lewej', a która 'po prawej'
                                            if (dr[0]>0) bAngle-=C;
                                            else aAngle-=C;
                                        }
                                        if (checkMoleculeInteriorsOverlap(drRoot,aAngle,bAngle)<0) collidingPairs++;
                                    }
                                }
                            }

                            rho=N/volume; pacFrac=1.0/VcpPerParticle/rho;
                            if (cycle%intervalResults==0)
                                fprintf(fileAllResults,"%.17E\t%.17E\t%.17E\t\n",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1]);

                            if (cycle%intervalOrientations==0) {
                                /////skan po konkretnej cząstce (np. dla 224: 14*(16/2)-(14/2)=105, etc.) - w srodku by nie skakala na granicy pudla periodycznego; bezposrednie uzycie w Mathematice (format tablicy)
                                if (cycle-cyclesOfEquilibrationBuffer>=cyclesOfMeasurementBuffer)
                                    fprintf(fileOrientations,"{%.12E,%.12E,%.12E}}",particles[indexScanned].r[0],particles[indexScanned].r[1],particles[indexScanned].phi);
                                else fprintf(fileOrientations,"{%.12E,%.12E,%.12E},",particles[indexScanned].r[0],particles[indexScanned].r[1],particles[indexScanned].phi);
                                /////skan po wszystkich cząstkach
                                for (int i=0;i<activeN-1;i++) fprintf(fileAllOrientations,"%.12E,",particles[i].phi);
                                fprintf(fileAllOrientations,"%.12E\n",particles[activeN-1].phi);
                                /////OCF
                                if (OCFMode) {
                                    char bufferText[4096]="", bufferAngle[20];
                                    for (int i=0;i<activeN;i++) for (int j=0;j<particles[i].neighCounter;j++) {
                                        if (i<particles[i].neighbours[j]) {
                                            strcpy(bufferText,"");
                                            getParticlesDistanceSquared(&particles[i],&particles[particles[i].neighbours[j]],boxMatrix);
                                            double gamma=atan(dr[1]/dr[0]),
                                                   aAngle=particles[particles[i].neighbours[j]].phi-gamma,
                                                   bAngle=particles[i].phi-gamma;
                                            if (multimerN%2!=0) {
                                                if (dr[0]>0) bAngle-=C;
                                                else aAngle-=C;
                                            }
                                            aAngle=normalizeAngle(aAngle+C); bAngle=normalizeAngle(bAngle+C);

                                            strncat(bufferText,"{",2); sprintf(bufferAngle,"%.12E",aAngle); strncat(bufferText,bufferAngle,20);
                                            strncat(bufferText,",",2); sprintf(bufferAngle,"%.12E",bAngle); strncat(bufferText,bufferAngle,20);
                                            if (i>=activeN-2 && cycle-cyclesOfEquilibrationBuffer>=cyclesOfMeasurementBuffer) strncat(bufferText,"}}",3); // -2 because it's last particle which has neighbour UNtested (the last one)
                                            else strncat(bufferText,"},",3);
                                            fprintf(fileOrientatCorrelFun,"%s",bufferText);
                                        }
                                    }
                                }
                            }

                            if (saveConfigurations && cycle%savedConfigurationsInt==0) {
                                fprintf(fileSavedConfigurations,"%ld\t%.12E\t%.12E\t%.12E\t{",(cycle+(long)args[4]),boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1]);
                                int activeNMinus1=activeN-1;
                                for (int i=0;i<activeNMinus1;i++)
                                    fprintf(fileSavedConfigurations,"m[%.17E,%.17E,%.17E],",particles[i].r[0],particles[i].r[1],particles[i].phi);
                                fprintf(fileSavedConfigurations,"m[%.17E,%.17E,%.17E]}",particles[activeNMinus1].r[0],particles[activeNMinus1].r[1],particles[activeNMinus1].phi);
                                fprintf(fileSavedConfigurations,"\n");
                            }

                            if (cycle%intervalOutput==0) {
                                if (countCollidingPairs) printf("Cycle: %ld, CollPairs: %d, ",(cycle+(long)args[4]),collidingPairs);
                                else printf("Cycle: %ld, ",(cycle+(long)args[4]));
                                printf("simulation time: full-%ldsec, measurement-%ldsec\n",time(0)-timeStart,time(0)-timeEquilibration);
                                printf("   AccRatR: %.4E, dR: %.4E, AccRatV: %.4E, dV: %.4E\n",acceptanceRatioR,deltaR,acceptanceRatioV,deltaV);
                                printf("   Dens: %.4E, V/V_cp: %.4E, PressRed: %.4E\n",rho,pacFrac,pressureReduced);
                                printf("   box00: %.8E, box11: %.8E, box01(10): %.8E\n",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1]);
                                printf("   totalInterEnergy: %.8E\n",totalInteractionEnergy);
                                fflush(stdout);

                                fileConfigurations = fopen(bufferConfigurations,"w");
                                fprintf(fileConfigurations,"Rho: %.12E\tV/V_cp: %.12E\tPressureRed: %.12E\tRandStart: %.1f\tRandSteps: %.1f\tCycles: %ld\tEquilTime: %ldsec\tMeasuTime: %ldsec\n\n",rho,pacFrac,
                                        pressureReduced,randomStartStep[0],randomStartStep[1],(cycle+(long)args[4]),(timeEquilibration-timeStart),(timeEnd-timeEquilibration));
                                fprintf(fileConfigurations,"=====================\tConfigurations (data to load)\t=====================\n");
                                fprintf(fileConfigurations,"%.1f %.1f %.17E %.17E %ld %.17E %.17E %.17E %.17E %.17E {",randomStartStep[0],randomStartStep[1],rho,pressureReduced,(cycle+(long)args[4]),boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1],deltaR,deltaV);
                                for (int i=0;i<activeN;i++)
                                    fprintf(fileConfigurations,"m[%.17E,%.17E,%.17E,%.2E,%.6E,%d],",particles[i].r[0],particles[i].r[1],particles[i].phi,multimerS,multimerD,multimerN);
                                fprintf(fileConfigurations,"{Opacity[0.2],Red,Polygon[{{0,0},{%.12E,%.12E},{%.12E,%.12E},{%.12E,%.12E}}]},{Opacity[0.2],Green,Disk[{%.12E,%.12E},%.12E]}}",boxMatrix[0][0],boxMatrix[1][0],boxMatrix[0][0]+boxMatrix[0][1],boxMatrix[1][0]+boxMatrix[1][1],boxMatrix[0][1],boxMatrix[1][1],particles[indexScanned].r[0],particles[indexScanned].r[1],neighRadius);
                                fprintf(fileConfigurations,"\n==========================================\n\n\nboxMatrix[0][0]=%.12E, boxMatrix[1][1]=%.12E, boxMatrix[0][1]=boxMatrix[1][0]=%.12E",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1]);
                                fclose(fileConfigurations);
                            }
                        } else {//dostosowywanie delt w trakcie pomiaru
                            if (acceptanceRatioR>desiredAcceptanceRatioR) deltaR*=1.05; else deltaR*=0.95;
                            if (deltaR>maxDeltaR*multimerD) deltaR=maxDeltaR*multimerD;  //jednostka odległości to d^**\sigma
                            deltaPhi=deltaR*2.0*sin(C)/multimerS;
                            if (acceptanceRatioV>desiredAcceptanceRatioV) deltaV*=1.05; else deltaV*=0.95;

                            int sampleNumberMod100=(cycle/intervalSampling)%100;
                            updateTableAndGetActualMean(deltaRTable,deltaRMean,sampleNumberMod100,deltaR); deltaR=deltaRMean;
                            updateTableAndGetActualMean(deltaVTable,deltaVMean,sampleNumberMod100,deltaV); deltaV=deltaVMean;
                        }
                            //printf("c: %ld, v*: %.17E, rho: %.17E\n\trR: %.17E -d/aR: %ld/%ld, dR: %.17E\n\trV: %.17E -d/aV: %ld/%ld, dV: %.17E\n\tb00: %.17E, b11: %.17E, b01: %.17E\n",cycle,volume/VcpPerParticle/N,N/volume,acceptanceRatioR,displacedNumberR,attemptedNumberR,deltaR,acceptanceRatioV,displacedNumberV,attemptedNumberV,deltaV,boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1]);
                        attemptedNumberR=0; displacedNumberR=0;
                        attemptedNumberV=0; displacedNumberV=0;
                    }
                    if (simulationStage==1 && cycle-cyclesOfEquilibrationBuffer>=cyclesOfMeasurementBuffer) simulationStage=2;
                }
            }
            fclose(fileAllResults);
            fclose(fileOrientations); fclose(fileAllOrientations);
            if (saveConfigurations) fclose(fileSavedConfigurations); if (OCFMode) fclose(fileOrientatCorrelFun);
            if (timeEquilibration==0) timeEquilibration=time(0);
            printf("Checking molecule's core overlaps in final configuration... "); fflush(stdout);
            if (!onlyMath[0] && getEnergyAll(particles,boxMatrix,1)<0) {
                printf("Final configuration contains molecule's core overlap(s) [energy->infinity]. ");
                char allResultsErrorFileName[200],configurationErrorFileName[200];
                strcpy(allResultsErrorFileName,allResultsFileName); strncat(allResultsErrorFileName,".err",5);
                strcpy(configurationErrorFileName,bufferConfigurations); strncat(configurationErrorFileName,".err",5);
                if (rename(allResultsFileName,allResultsErrorFileName)==0) printf("Results file successfully renamed (.err). ");
                else printf("Error renaming results file (.err). ");
                if (rename(bufferConfigurations,configurationErrorFileName)==0) printf("Configuration file successfully renamed (.err).\n");
                else printf("Error renaming configuration file (.err).\n");
                return 0;
            } else {printf("done\n"); fflush(stdout);}
            timeEnd=time(0);




/////////////////////////////////////////////// OBLICZENIE WYNIKOW

            printf("Start of calculation of results...\n");

            //obliczenie liczby linii danych (potrzebne do podziału na zespoły i obliczenia średnich błędów)
            printf("Calculation of data lines... "); fflush(stdout);
            fileAllResults=fopen(allResultsFileName,"rt");
            char linia[4096]; double dataLicznik=0; int faultyLines=0, onlyMathLinesBuffer=onlyMath[1];
            if (onlyMath[0]) for (int i=0;i<onlyMathLinesBuffer;i++) {
                if (fgets(linia,300,fileAllResults)!=NULL && !isLineCorrect(linia)) onlyMathLinesBuffer++;
            }
            while(fgets(linia,300,fileAllResults)!=NULL) {
                if (isLineCorrect(linia)) dataLicznik++; else faultyLines++;
            }
            fclose(fileAllResults);
            printf("done (Found %ld data lines [%d faulty lines occurred].",(long)dataLicznik,faultyLines);
            if ((long)dataLicznik%10>0) printf(" Last %ld lines won't be considered, due to calculations of averages in 10 sets.)\n",(long)dataLicznik%10); else printf(")\n");
            dataLicznik-=(long)dataLicznik%10;

            //obliczenie srednich wartosci mierzonych wielkosci
            printf("Calculation of averages... "); fflush(stdout);
            double avVolumeSet[10], avBoxMatrixSet[10][3], avRhoSet[10], avPacFracSet[10]; //wyniki dzielone są na 10 zespołów (obliczane są nieskorelowane wzajemnie średnie "lokalne", do obliczenia błędu średniej "globalnej")
            for (int i=0;i<10;i++) {
                avVolumeSet[i]=0; for (int j=0;j<3;j++) avBoxMatrixSet[i][j]=0;
                avRhoSet[i]=0; avPacFracSet[i]=0;
            }
            fileAllResults=fopen(allResultsFileName,"rt");
            if (onlyMath[0]) for (int i=0;i<onlyMathLinesBuffer;i++) fgets(linia,300,fileAllResults);
            double lineCounter=0;
            while(fgets(linia,300,fileAllResults)!=NULL && lineCounter<dataLicznik) {
                sscanf(linia,"%c",linia);
                int actIndex=0;
                int dataIndex=0; double dataD[3]; while (dataIndex<3) {
                    char data[50]="";
                    int licznik=0, dotCounter=0;
                    while (linia[actIndex]!='\t' && (dotCounter=linia[actIndex]=='.'?dotCounter+1:dotCounter)<=1 && licznik<50) data[licznik++]=linia[actIndex++];
                    if (dotCounter>1 || licznik>=50) {dataIndex=10; continue;}
                    actIndex++;
                    if (dataIndex<10 && ((data[0]!='-' && data[1]!='.') || (data[0]=='-' && data[2]!='.'))) dataIndex=10;
                    else dataD[dataIndex++]=strtod(data,NULL);
                } if (dataIndex<10 && linia[actIndex]!='\n') dataIndex=10;
                if (dataIndex<10) {
                    int setIndex=(int)(lineCounter/dataLicznik*10.0);
                    for (int i=0;i<3;i++) avBoxMatrixSet[setIndex][i]+=dataD[i];
                    lineCounter++;
                }
            }
            fclose(fileAllResults);
            double avVolume=0, avBoxMatrix[3]={0,0,0}, avRho=0, avPacFrac=0;
            for (int i=0;i<10;i++) {
                for (int j=0;j<3;j++) {
                    avBoxMatrixSet[i][j]/=dataLicznik*0.1; avBoxMatrix[j]+=avBoxMatrixSet[i][j];
                }
                avVolumeSet[i]=fabs(avBoxMatrixSet[i][0]*avBoxMatrixSet[i][1]-avBoxMatrixSet[i][2]*avBoxMatrixSet[i][2]); avVolume+=avVolumeSet[i];
                avRhoSet[i]=N/avVolumeSet[i]; avRho+=avRhoSet[i];
                avPacFracSet[i]=1.0/VcpPerParticle/avRhoSet[i]; avPacFrac+=avPacFracSet[i];
            }
            avVolume*=0.1; avRho*=0.1; avPacFrac*=0.1; for (int i=0;i<3;i++) avBoxMatrix[i]*=0.1;
            //obliczenie bledow mierzonych wielkosci
            double dAvVolume=0, dAvBoxMatrix[3]={0,0,0}, dAvRho=0, dAvPacFrac=0;
            for (int i=0;i<10;i++) {double epsilon=avVolume-avVolumeSet[i]; dAvVolume+=epsilon*epsilon;} dAvVolume=getAvErrorFromSumEps(dAvVolume,90.0); //10*9 (n(n-1))
            for (int j=0;j<3;j++) {for (int i=0;i<10;i++) {double epsilon=avBoxMatrix[j]-avBoxMatrixSet[i][j]; dAvBoxMatrix[j]+=epsilon*epsilon;} dAvBoxMatrix[j]=getAvErrorFromSumEps(dAvBoxMatrix[j],90.0);}
            for (int i=0;i<10;i++) {double epsilon=avRho-avRhoSet[i]; dAvRho+=epsilon*epsilon;} dAvRho=getAvErrorFromSumEps(dAvRho,90.0);
            for (int i=0;i<10;i++) {double epsilon=avPacFrac-avPacFracSet[i]; dAvPacFrac+=epsilon*epsilon;} dAvPacFrac=getAvErrorFromSumEps(dAvPacFrac,90.0);
            printf("done\n");

            //obliczenie srednich iloczynow elementow tensora odkształceń
            printf("Calculation of average values of products of strain tensor's elements... "); fflush(stdout);
            double e1111Set[10], e1122Set[10], e1212Set[10], e2222Set[10], e1112Set[10], e1222Set[10],
                   HxyHyx,HxxHyy,HxxHxy,Hxx2,Hyy2,mod0,mod1;
            for (int i=0;i<10;i++) {e1111Set[i]=0; e1122Set[i]=0; e1212Set[i]=0; e2222Set[i]=0; e1112Set[i]=0; e1222Set[i]=0;}
            fileAllResults=fopen(allResultsFileName,"rt");
            if (onlyMath[0]) for (int i=0;i<onlyMathLinesBuffer;i++) fgets(linia,300,fileAllResults);
            lineCounter=0; int setIndex, oldSetIndex=-1;
            while(fgets(linia,300,fileAllResults)!=NULL && lineCounter<dataLicznik) {
                setIndex=(int)(lineCounter/dataLicznik*10.0);
                if (setIndex!=oldSetIndex) {
                    HxyHyx=avBoxMatrixSet[setIndex][2]*avBoxMatrixSet[setIndex][2];
                    HxxHyy=avBoxMatrixSet[setIndex][0]*avBoxMatrixSet[setIndex][1];
                    HxxHxy=avBoxMatrixSet[setIndex][0]*avBoxMatrixSet[setIndex][2];
                    Hxx2=avBoxMatrixSet[setIndex][0]*avBoxMatrixSet[setIndex][0];
                    Hyy2=avBoxMatrixSet[setIndex][1]*avBoxMatrixSet[setIndex][1];
                    mod0=HxyHyx-HxxHyy;
                    mod1=1.0/(2.0*mod0*mod0);
                    oldSetIndex=setIndex;
                }
                sscanf(linia,"%c",linia);
                int actIndex=0;
                double h[3],h11,h22,h12;
                int dataIndex=0; while (dataIndex<3) {
                    char data[50]="";
                    int licznik=0, dotCounter=0;
                    while (linia[actIndex]!='\t' && (dotCounter=linia[actIndex]=='.'?dotCounter+1:dotCounter)<=1 && licznik<50) data[licznik++]=linia[actIndex++];
                    if (dotCounter>1 || licznik>=50) {dataIndex=10; continue;}
                    actIndex++;
                    if (dataIndex<10) {
                        if ((data[0]!='-' && data[1]!='.') || (data[0]=='-' && data[2]!='.')) dataIndex=10;
                        else h[dataIndex++]=strtod(data,NULL);
                    }
                } if (dataIndex<10 && linia[actIndex]!='\n') dataIndex=10;
                if (dataIndex<10) {
                    h11=h[0]; h22=h[1]; h12=h[2];
                    double hxyhyx=h12*h12,
                           hxxPhyy=h11+h22,
                           hxx2=h11*h11,
                           hyy2=h22*h22,

                           e11=mod1*(HxyHyx*(hxyhyx-HxyHyx+hyy2)-2.0*(h11*avBoxMatrixSet[setIndex][2]*h12-avBoxMatrixSet[setIndex][0]*HxyHyx+avBoxMatrixSet[setIndex][2]*h12*h22)*avBoxMatrixSet[setIndex][1]+(hxx2-Hxx2+hxyhyx)*Hyy2),
                           e22=mod1*(HxyHyx*(hxx2+hxyhyx-HxyHyx)-2.0*(HxxHxy*h12*hxxPhyy-HxxHyy*HxyHyx)+Hxx2*(hxyhyx+hyy2-Hyy2)),
                           e12=mod1*(-HxxHxy*(hxyhyx+hyy2)+HxxHyy*h12*hxxPhyy+avBoxMatrixSet[setIndex][2]*(h12*avBoxMatrixSet[setIndex][2]*hxxPhyy-(hxx2+hxyhyx)*avBoxMatrixSet[setIndex][1]));

                    e1111Set[setIndex]+=e11*e11;
                    e1122Set[setIndex]+=e11*e22;
                    e1212Set[setIndex]+=e12*e12;
                    e2222Set[setIndex]+=e22*e22;
                    e1112Set[setIndex]+=e11*e12;
                    e1222Set[setIndex]+=e12*e22;
                    lineCounter++;
                }
            }
            fclose(fileAllResults);
            double e1111=0, e1122=0, e1212=0, e2222=0, e1112=0, e1222=0;
            for (int i=0;i<10;i++) {
                e1111Set[i]/=dataLicznik*0.1; e1111+=e1111Set[i];
                e1122Set[i]/=dataLicznik*0.1; e1122+=e1122Set[i];
                e1212Set[i]/=dataLicznik*0.1; e1212+=e1212Set[i];
                e2222Set[i]/=dataLicznik*0.1; e2222+=e2222Set[i];
                e1112Set[i]/=dataLicznik*0.1; e1112+=e1112Set[i];
                e1222Set[i]/=dataLicznik*0.1; e1222+=e1222Set[i];
            }
            e1111*=0.1; e1122*=0.1; e1212*=0.1; e2222*=0.1; e1112*=0.1; e1222*=0.1;
            //obliczenie bledow iloczynow elementow tensora odkształceń
            double dE1111=0, dE1122=0, dE1212=0, dE2222=0, dE1112=0, dE1222=0;
            for (int i=0;i<10;i++) {double epsilon=e1111-e1111Set[i]; dE1111+=epsilon*epsilon;} dE1111=getAvErrorFromSumEps(dE1111,90.0); //10*9 (n(n-1))
            for (int i=0;i<10;i++) {double epsilon=e1122-e1122Set[i]; dE1122+=epsilon*epsilon;} dE1122=getAvErrorFromSumEps(dE1122,90.0);
            for (int i=0;i<10;i++) {double epsilon=e1212-e1212Set[i]; dE1212+=epsilon*epsilon;} dE1212=getAvErrorFromSumEps(dE1212,90.0);
            for (int i=0;i<10;i++) {double epsilon=e2222-e2222Set[i]; dE2222+=epsilon*epsilon;} dE2222=getAvErrorFromSumEps(dE2222,90.0);
            for (int i=0;i<10;i++) {double epsilon=e1112-e1112Set[i]; dE1112+=epsilon*epsilon;} dE1112=getAvErrorFromSumEps(dE1112,90.0);
            for (int i=0;i<10;i++) {double epsilon=e1222-e1222Set[i]; dE1222+=epsilon*epsilon;} dE1222=getAvErrorFromSumEps(dE1222,90.0);
            printf("done\n");

            //obliczenie podatnosci, wspolczynnika Poissona i modulow sprezystosci
            printf("Calculation of compliances, Poisson's ratio and elastic moduli... "); fflush(stdout);
            //Eijkl - bezwymiarowe [strain], volume - [sigma^2], T*=kT/eps [eps - jednostka energii (w domyśle =1)]
            double s1111=e1111*avVolume/TReduced, dS1111=(fabs(e1111*dAvVolume)+fabs(dE1111*avVolume))/TReduced,
                   s1122=e1122*avVolume/TReduced, dS1122=(fabs(e1122*dAvVolume)+fabs(dE1122*avVolume))/TReduced,
                   s1212=e1212*avVolume/TReduced, dS1212=(fabs(e1212*dAvVolume)+fabs(dE1212*avVolume))/TReduced,
                   s2222=e2222*avVolume/TReduced, dS2222=(fabs(e2222*dAvVolume)+fabs(dE2222*avVolume))/TReduced,
                   s1112=e1112*avVolume/TReduced, dS1112=(fabs(e1112*dAvVolume)+fabs(dE1112*avVolume))/TReduced,
                   s1222=e1222*avVolume/TReduced, dS1222=(fabs(e1222*dAvVolume)+fabs(dE1222*avVolume))/TReduced,

                   S11=(s1111+s2222)*0.5,
                   S66=4.0*s1212,
                   avNu=-s1122/S11, dAvNu=fabs(dS1122/S11)+fabs((dS1111+dS2222)*0.5*s1122/S11/S11), //nu obliczane z Sxxyy
                   avNu2=S66/S11*0.5-1, dAvNu2=fabs(0.5/S11*4.0*dS1212)+fabs(0.5*S66/S11/S11*(dS1111+dS2222)*0.5), //nu obliczane z Sxyxy (inny rodzaj scinania, przy izotropowych ukladach powinno byc tyle samo co nu1)

                   //stary sposob obliczania nu (nu_x, nu_y i srednia)
                   /*nu2211_1111=-s1122/s1111, dNu2211_1111=fabs(dS1122/s1111)+fabs(dS1111*s1122/s1111/s1111),
                   nu1122_2222=-s1122/s2222, dNu1122_2222=fabs(dS1122/s2222)+fabs(dS2222*s1122/s2222/s2222),
                   avNu=(nu2211_1111+nu1122_2222)/2.0, dAvNu=(dNu2211_1111+dNu1122_2222)/2.0,*/

                   //\lambdaReduced=\lambda*d^**\sigma^2/kT, T*=kT/eps [eps - jednostka energii (w domyśle =1)]
                   l11=multimerD*multimerS/TReduced/(8.0*(s1111+s1122)), dL11=multimerD*multimerS/TReduced*(fabs(dS1111)+fabs(dS1122))/(8.0*fabs(s1111+s1122)*fabs(s1111+s1122)),
                   l12=multimerD*multimerS/TReduced/(8.0*(s2222+s1122)), dL12=multimerD*multimerS/TReduced*(fabs(dS2222)+fabs(dS1122))/(8.0*fabs(s2222+s1122)*fabs(s2222+s1122)),

                   B1=4.0*l11, dB1=4.0*dL11,
                   B2=4.0*l12, dB2=4.0*dL12,
                   avB=(B1+B2)/2.0, dAvB=(dB1+dB2)/2.0,

                   my1=(B1-B1*avNu)/(1.0+avNu), dMy1=fabs(dB1*(1.0-avNu)/(1.0+avNu))+fabs(dAvNu*(-B1/(1.0+avNu)-(B1-B1*avNu)/(1.0+avNu)/(1.0+avNu))),
                   my2=(B2-B2*avNu)/(1.0+avNu), dMy2=fabs(dB2*(1.0-avNu)/(1.0+avNu))+fabs(dAvNu*(-B2/(1.0+avNu)-(B2-B2*avNu)/(1.0+avNu)/(1.0+avNu))),
                   avMy=(my1+my2)/2.0, dAvMy=(dMy1+dMy2)/2.0,

                   avE=4.0*avB*avMy/(avB+avMy), dAvE=4.0*(fabs(avB*avB*dAvMy)+fabs(dAvB*avMy*avMy))/fabs(avB+avMy)/fabs(avB+avMy);
            printf("done\n");

            //tworzenie pliku wynikowego do Origina - rozklad orientacyjny dla 1 czastki
            printf("Creation of 1-particle orientation file for Origin... "); fflush(stdout);
            double componentCounter=0, averageCos6PhiOne=0, ODFMaxOne=0;
            fileOrientations=fopen(bufferOrientations,"rt");
            double ODF_1P[ODFLength]; for (int i=0;i<ODFLength;i++) ODF_1P[i]=0; //Orientational Distribution Function (1 Particle)
            int licznik=1,character;
            while ((character=fgetc(fileOrientations))!=EOF) {
                if (character==',') licznik++;
                if (licznik==3) {
                    licznik=0;
                    char data[50]=""; int actIndex=0;
                    while (true) {
                        character=fgetc(fileOrientations);
                        if (character!='}' && actIndex<50) data[actIndex++]=character;
                        else {
                            double angle = normalizeAngle(strtod(data,NULL)+C);
                            averageCos6PhiOne+=cos(6.0*angle); componentCounter++;
                            int index = round((angle+C)/2.0/C*(double)(ODFLength-1.0));
                            ODF_1P[index]++;
                            break;
                        }
                    }
                }
            }
            fclose(fileOrientations);
            double dPhi=2.0*C/((double)(ODFLength-1.0));
            double suma=0; for (int i=0;i<ODFLength;i++) suma+=ODF_1P[i]; for (int i=0;i<ODFLength;i++) ODF_1P[i]/=suma*dPhi;
            averageCos6PhiOne/=componentCounter; for (int i=0;i<ODFLength;i++) if (ODFMaxOne<ODF_1P[i]) ODFMaxOne=ODF_1P[i];
            printf("done\n");

            //tworzenie pliku wynikowego do Origina - rozklad orientacyjny dla wszystkich czastek
            printf("Creation of ALL-particle orientation file for Origin... "); fflush(stdout);
            componentCounter=0; double averageCos6PhiAll=0, ODFMaxAll=0, averagePhiAll=0;
            fileAllOrientations = fopen(allOrientationsFileName,"rt");
            double ODF_AllP[ODFLength]; for (int i=0;i<ODFLength;i++) ODF_AllP[i]=0; //Orientational Distribution Function (All Particles)
            char data[50]=""; int actIndex=0;
            while ((character=fgetc(fileAllOrientations))!=EOF) {
                if (character!=',' && character!='\n') data[actIndex++]=character;
                else {
                    data[actIndex++]=' '; actIndex=0;
                    double angle = normalizeAngle(strtod(data,NULL)+C);
                    averagePhiAll+=angle;
                    averageCos6PhiAll+=cos(6.0*angle); componentCounter++;
                    int index = round((angle+C)/2.0/C*(double)(ODFLength-1.0));
                    ODF_AllP[index]++;
                }
            }
            fclose(fileAllOrientations);
            suma=0; for (int i=0;i<ODFLength;i++) suma+=ODF_AllP[i]; for (int i=0;i<ODFLength;i++) ODF_AllP[i]/=suma*dPhi;
            averagePhiAll/=componentCounter; averageCos6PhiAll/=componentCounter;
            int maxODFAllIndex;
            for (int i=0;i<ODFLength;i++) if (ODFMaxAll<ODF_AllP[i]) {
                ODFMaxAll=ODF_AllP[i];
                maxODFAllIndex=i;
            }
            printf("done\n");

            //analiza konfiguracji przejsciowych (dPhi z konfiguracji na konfiguracje)
            double avAbsDPhi=0;
            if (saveConfigurations) {
                printf("Transient configurations analysis... "); fflush(stdout);
                fileSavedConfigurations = fopen(bufferSavedConfigurations,"rt");
                double prevCfg[activeN][3]; bool undefinedPrevCfg=true;
                int CfgQuantity=0,CfgFileQuantity=0;
                while ((character=fgetc(fileSavedConfigurations))!=EOF) {
                    if (character=='n') {//text: newCfgFile (after merging) [it's NOT \n -> it's n, first char of 'new']
                        undefinedPrevCfg=true;
                        while (fgetc(fileSavedConfigurations)!='\n');
                        continue;
                    } else CfgQuantity++;
                    while (fgetc(fileSavedConfigurations)!='[');
                    for (int i=0;i<activeN;i++) {
                        for (int j=0;j<3;j++) {
                            strcpy(data,""); actIndex=0;
                            while ((character=fgetc(fileSavedConfigurations))!=',' && character!=']') data[actIndex++]=character;
                            data[actIndex++]=' ';
                            if (j==2) {
                                for (int k=0;k<3;k++) fgetc(fileSavedConfigurations); //skip ",m["
                                if (!undefinedPrevCfg) {
                                    avAbsDPhi+=fabs(strtod(data,NULL)-prevCfg[i][j]);
                                }
                            }
                            prevCfg[i][j]=strtod(data,NULL);
                        }
                    }
                    if (undefinedPrevCfg) {undefinedPrevCfg=false; CfgFileQuantity++;}
                }
                fclose(fileSavedConfigurations);
                avAbsDPhi/=(double)activeN*(CfgQuantity-CfgFileQuantity);
                printf("done\n");
            }

            long timeEndMath=time(0);




/////////////////////////////////////////////// ZAPIS DANYCH DO PLIKU

            printf("Saving data to files... "); fflush(stdout);
            timeEq+=(timeEquilibration-timeStart); timeMe+=(timeEnd-timeEquilibration); timeMath+=(timeEndMath-timeEnd);

            fileResults = fopen(resultsFileName,"a");
            fileExcelResults = fopen(excelResultsFileName,"a");
            if (!onlyMath[0]) {
                fileConfigurations = fopen(bufferConfigurations,"w");
                fileConfigurationsList = fopen(configurationsListFileName,"a");
            }
            fileOrientationsResults = fopen(bufferOrientationsResults,"w");
            fileAllOrientationsResults = fopen(allOrientationsResultsFileName,"w");

            if (saveConfigurations) {
                fprintf(fileResults,"%ld\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%ld\t%.17E\n",(cycle+(long)args[4]),pressureReduced,avVolume,dAvVolume,avBoxMatrix[0],dAvBoxMatrix[0],avBoxMatrix[1],dAvBoxMatrix[1],avBoxMatrix[2],dAvBoxMatrix[2],avRho,dAvRho,avPacFrac,dAvPacFrac,s1111,dS1111,s1122,dS1122,s1212,dS1212,s2222,dS2222,s1112,dS1112,s1222,dS1222,avNu,dAvNu,avNu2,dAvNu2,avB,dAvB,avMy,dAvMy,avE,dAvE,ODFMaxOne,averageCos6PhiOne,ODFMaxAll,averageCos6PhiAll,-C+maxODFAllIndex*dPhi,averagePhiAll,savedConfigurationsInt,avAbsDPhi);
                fprintf(fileExcelResults,"%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%ld\t%.17E\n",pressureReduced,avPacFrac,avNu,avNu2,ODFMaxAll,averageCos6PhiAll,-C+maxODFAllIndex*dPhi,averagePhiAll,avB,avMy,avE,savedConfigurationsInt,avAbsDPhi);
            } else {
                fprintf(fileResults,"%ld\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\n",(cycle+(long)args[4]),pressureReduced,avVolume,dAvVolume,avBoxMatrix[0],dAvBoxMatrix[0],avBoxMatrix[1],dAvBoxMatrix[1],avBoxMatrix[2],dAvBoxMatrix[2],avRho,dAvRho,avPacFrac,dAvPacFrac,s1111,dS1111,s1122,dS1122,s1212,dS1212,s2222,dS2222,s1112,dS1112,s1222,dS1222,avNu,dAvNu,avNu2,dAvNu2,avB,dAvB,avMy,dAvMy,avE,dAvE,ODFMaxOne,averageCos6PhiOne,ODFMaxAll,averageCos6PhiAll,-C+maxODFAllIndex*dPhi,averagePhiAll);
                fprintf(fileExcelResults,"%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\n",pressureReduced,avPacFrac,avNu,avNu2,ODFMaxAll,averageCos6PhiAll,-C+maxODFAllIndex*dPhi,averagePhiAll,avB,avMy,avE);
            }

            if (!onlyMath[0]) {
                rho=N/volume; pacFrac=1.0/VcpPerParticle/rho;
                fprintf(fileConfigurations,"Rho: %.12E\tV/V_cp: %.12E\tPressureRed: %.12E\tRandStart: %.1f\tRandSteps: %.1f\tCycles: %ld\tEquilTime: %ldsec\tMeasuTime: %ldsec\n\n",rho,pacFrac,
                        pressureReduced,randomStartStep[0],randomStartStep[1],(cycle+(long)args[4]),(timeEquilibration-timeStart),(timeEnd-timeEquilibration));
                fprintf(fileConfigurations,"=====================\tConfigurations (data to load)\t=====================\n");
                fprintf(fileConfigurations,"%.1f %.1f %.17E %.17E %ld %.17E %.17E %.17E %.17E %.17E {",randomStartStep[0],randomStartStep[1],rho,pressureReduced,(cycle+(long)args[4]),boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1],deltaR,deltaV);
                fprintf(fileConfigurationsList,"multimers[x_,y_,kI_]:={");
                for (int i=0;i<activeN;i++) {
                    fprintf(fileConfigurations,"m[%.17E,%.17E,%.17E,%.2E,%.6E,%d],",particles[i].r[0],particles[i].r[1],particles[i].phi,multimerS,multimerD,multimerN);
                    fprintf(fileConfigurationsList,"m[%.12E+x,%.12E+y,%.12E,%.2E,%.6E,%d],",particles[i].r[0],particles[i].r[1],particles[i].phi,multimerS,multimerD,multimerN);
                }
                fprintf(fileConfigurations,"{Opacity[0.2],Red,Polygon[{{0,0},{%.12E,%.12E},{%.12E,%.12E},{%.12E,%.12E}}]},{Opacity[0.2],Green,Disk[{%.12E,%.12E},%.12E]}}",boxMatrix[0][0],boxMatrix[1][0],boxMatrix[0][0]+boxMatrix[0][1],boxMatrix[1][0]+boxMatrix[1][1],boxMatrix[0][1],boxMatrix[1][1],particles[indexScanned].r[0],particles[indexScanned].r[1],neighRadius);
                fprintf(fileConfigurations,"\n==========================================\n\n\nboxMatrix[0][0]=%.12E, boxMatrix[1][1]=%.12E, boxMatrix[0][1]=boxMatrix[1][0]=%.12E",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1]);
                fprintf(fileConfigurationsList,"{Opacity[If[x==0 && y==0,0.4,0]],Red,Polygon[{{0,0},{%.12E,%.12E},{%.12E,%.12E},{%.12E,%.12E}}]},{Opacity[If[x==0 && y==0,0.4,0]],Green,Disk[{%.12E,%.12E},%.12E]}};\nconfigurationsList=Append[configurationsList,g[%.12E,multimers[%.12E,%.12E,kolorIndex=1],multimers[%.12E,%.12E,kolorIndex=1],multimers[%.12E,%.12E,kolorIndex=1],multimers[0,0,kolorIndex=1]]];\n",
                        boxMatrix[0][0],boxMatrix[1][0],boxMatrix[0][0]+boxMatrix[0][1],boxMatrix[1][0]+boxMatrix[1][1],boxMatrix[0][1],boxMatrix[1][1],particles[indexScanned].r[0],particles[indexScanned].r[1],neighRadius,pacFrac,boxMatrix[0][0],boxMatrix[1][0],boxMatrix[0][0]+boxMatrix[0][1],boxMatrix[1][0]+boxMatrix[1][1],boxMatrix[0][1],boxMatrix[1][1]);
            }

            for (int i=0;i<ODFLength;i++) {
                fprintf(fileOrientationsResults,"%.12E\t%.12E\n",-C+i*dPhi,ODF_1P[i]);
                fprintf(fileAllOrientationsResults,"%.12E\t%.12E\n",-C+i*dPhi,ODF_AllP[i]);
            }

            fclose(fileResults); fclose(fileExcelResults);
            if (!onlyMath[0]) {
                fclose(fileConfigurations); fclose(fileConfigurationsList);
            }
            fclose(fileOrientationsResults); fclose(fileAllOrientationsResults);
            printf("done\n\n");
        }




/////////////////////////////////////////////// PRZYGOTOWANIE DO KOLEJNEJ ITERACJI

        getNextArgument(arg,true);
        for (int i=0;i<2;i++) for (int j=0;j<2;j++) oldBoxMatrix[i][j]=boxMatrix[i][j];
        oldTotalInteractionEnergy=totalInteractionEnergy;
        args[4]=0;
        loadedConfiguration=0;
        generatorStartPoint=0;
    }
    printf("\nTime for equilibrations: %ldsec, time for measurments: %ldsec, time for math: %ldsec.\n",timeEq,timeMe,timeMath);
    printf("\nSimulation has been completed.\n"); fflush(stdout);
}
