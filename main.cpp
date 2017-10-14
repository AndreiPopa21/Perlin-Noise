#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>

using namespace std;
ifstream f("perlin.in");
ofstream g("perlin.out");
ofstream h("perl.out");
ofstream e("perl1.out");
ofstream k("perl2.out");
ofstream l("perl3.out");

struct Vect2D{float x;float y;};


///loc de imbunatatiri la algoritmul de tratare a chunckurilor
///introdu octave

struct Point  {
                  int x;
                  int y;

                  bool isLattice;

                  Vect2D gradient;

                  ///used only when isLattice is false
                  Vect2D p00; ///top-left corner
                  Vect2D p01; ///top-right corner
                  Vect2D p10; ///bottom-left corner
                  Vect2D p11; ///bottom-right corner

                  ///the displacement vectors  Gradient-this.Vect2D
                  Vect2D disp00;
                  Vect2D disp01;
                  Vect2D disp10;
                  Vect2D disp11;  ///used along the gradient vector of the lattice points

                  ///the Q-s for each corner Q= Disp*G which are going to be used for the interpolation
                  double q00;
                  double q01;
                  double q10;
                  double q11;

                  double finalHeight;

              };


int HashTable[201];

Point MatrixOfVertices[1025][1025];


int frequency;
int length;

void AssignHashTable()
{

  for(int i=-100;i<=100;i++)
    HashTable[i+100]=i;

/*for(int i=0;i<256;i++)
    cout<<HashTable[i]<<' ';*/

    for(int i=200;i>=0;i--) ///shuffle function
    {
      int t= rand()%(i+1);
      int r= HashTable[i];
      HashTable[i]=HashTable[t];
      HashTable[t]=r;
    }



}

void AssignPseudGradVec()
{
    int i,j;
    for(i=0;i<length+1;i++)
        {
            for(j=0;j<length+1;j++)
            {
                int t=i+j;
                int hashI= ((t*989)+89+i)%201;
                int hashJ= ((t*781)+56+j)%201;

                MatrixOfVertices[i][j].gradient.x=((HashTable[hashI]%100)*0.01+HashTable[hashI]/100);
                MatrixOfVertices[i][j].gradient.y=((HashTable[hashJ]%100)*0.01+HashTable[hashJ]/100); ///merge

               /* g<<"Gradientii punctului:  "<<i<<' '<<j<<"  sunt: "<<'\n';
                g<<"("<<MatrixOfVertices[i][j].gradient.x<<";"<<MatrixOfVertices[i][j].gradient.y<<")"<<'\n';
                g<<'\n';*/
            }


        }

    /*cout<<"("<<MatrixOfVertices[0][0].gradient.x<<";"<<MatrixOfVertices[0][0].gradient.y<<")"<<'\n';
    cout<<"("<<MatrixOfVertices[0][4].gradient.x<<";"<<MatrixOfVertices[0][4].gradient.y<<")"<<'\n';
    cout<<"("<<MatrixOfVertices[4][0].gradient.x<<";"<<MatrixOfVertices[4][0].gradient.y<<")"<<'\n';
    cout<<"("<<MatrixOfVertices[4][4].gradient.x<<";"<<MatrixOfVertices[4][4].gradient.y<<")"<<'\n';
*/
}

void AssignLatticePoints(int freq)
{
    int i,j;
    for(i=0;i<length+1;i=i+freq)
        for(j=0;j<length+1;j=j+freq)
        MatrixOfVertices[i][j].isLattice=true;

}

double NonLinearInterpolation(double t)
{

    return t * t * t * (t * (t * 6 - 15) + 10);
}

void DotProduct(int i,int j)
{
    Point t= MatrixOfVertices[i][j];
    ///p00
    t.q00=t.disp00.x*MatrixOfVertices[(int)t.p00.x][(int)t.p00.y].gradient.x+t.disp00.y*MatrixOfVertices[(int)t.p00.x][(int)t.p00.y].gradient.y;

    ///p01
    t.q01=t.disp01.x*MatrixOfVertices[(int)t.p01.x][(int)t.p01.y].gradient.x+t.disp01.y*MatrixOfVertices[(int)t.p01.x][(int)t.p01.y].gradient.y;

    ///p10
    t.q10=t.disp10.x*MatrixOfVertices[(int)t.p10.x][(int)t.p10.y].gradient.x+t.disp10.y*MatrixOfVertices[(int)t.p10.x][(int)t.p10.y].gradient.y;

    ///p11
    t.q11=t.disp11.x*MatrixOfVertices[(int)t.p11.x][(int)t.p11.y].gradient.x+t.disp11.y*MatrixOfVertices[(int)t.p11.x][(int)t.p11.y].gradient.y;

    MatrixOfVertices[i][j].q00=t.q00;
    MatrixOfVertices[i][j].q01=t.q01;
    MatrixOfVertices[i][j].q10=t.q10;
    MatrixOfVertices[i][j].q11=t.q11;

    /*k<<"Valorile finale ale punctului:  "<<i<<' '<<j<<"  sunt: "<<'\n';
    k<<MatrixOfVertices[i][j].q00<<' '<<MatrixOfVertices[i][j].q01<<' '<<MatrixOfVertices[i][j].q10<<' '<<MatrixOfVertices[i][j].q11<<endl;
    k<<'\n';*/

}

double lerp(double u,double v,double t)
{
    return u*(1-t)+v*t;
}

double interpPos(int x,int y,int t)
{
    double kk=((double)t-(double)x)/((double)y-(double)x);


  return kk;
}



double mapFunction(double value, double linf,double lsup, double dinf,double dsup)
{
    double tt=(value-linf)/(lsup-linf);
    return tt * tt * tt * (tt * (tt * 6 - 15) + 10);
}

void Fade(int x,int y)
{
   int pointTopRightX=MatrixOfVertices[x][y].p11.x;
   int pointTopRightY=MatrixOfVertices[x][y].p11.y;

   int pointTopLeftX=MatrixOfVertices[x][y].p10.x;
   int pointTopLeftY=MatrixOfVertices[x][y].p10.y;

   int pointBottomLeftX=MatrixOfVertices[x][y].p00.x;
   int pointBottomLeftY=MatrixOfVertices[x][y].p00.y;

   int pointBottomRightX=MatrixOfVertices[x][y].p01.x;
   int pointBottomRightY=MatrixOfVertices[x][y].p01.y;

   double tX1=interpPos(pointBottomLeftX,pointTopLeftX,x);
   double tX2=interpPos(pointBottomRightX,pointTopRightX,x);
   double tY=interpPos(pointTopLeftY,pointTopRightY,y);



   double s1=lerp(MatrixOfVertices[x][y].q00,MatrixOfVertices[x][y].q10,tX1);
   double s2=lerp(MatrixOfVertices[x][y].q01,MatrixOfVertices[x][y].q11,tX2);

   double r= NonLinearInterpolation(tY);


   MatrixOfVertices[x][y].finalHeight=lerp(s1,s2,r);

}

void PickSamplePoint(int freq)
{
    int i,j;
    for(i=0;i<length;i=i+freq)
        for(j=0;j<length;j=j+freq)
       {
        int s,t;
        for(s=i;s<i+freq;s++)
            for(t=j;t<j+freq;t++)
            {
                ///assign the the corner points of its lattice
                MatrixOfVertices[s][t].p00.x=i;
                MatrixOfVertices[s][t].p00.y=j;

                MatrixOfVertices[s][t].p01.x=i;
                MatrixOfVertices[s][t].p01.y=j+freq;

                MatrixOfVertices[s][t].p10.x=i+freq;
                MatrixOfVertices[s][t].p10.y=j;

                MatrixOfVertices[s][t].p11.x=i+freq;
                MatrixOfVertices[s][t].p11.y=j+freq;


                /*h<<MatrixOfVertices[s][t].p00.x<<'\t'<<MatrixOfVertices[s][t].p00.y<<'\n';
                h<<MatrixOfVertices[s][t].p01.x<<'\t'<<MatrixOfVertices[s][t].p01.y<<'\n';
                h<<MatrixOfVertices[s][t].p10.x<<'\t'<<MatrixOfVertices[s][t].p10.y<<'\n';
                h<<MatrixOfVertices[s][t].p11.x<<'\t'<<MatrixOfVertices[s][t].p11.y<<'\n';

                 h <<endl;*/

                ///assign the displacement of the corner points relative to the sample point

                MatrixOfVertices[s][t].disp00.x=s-MatrixOfVertices[s][t].p00.x;
                MatrixOfVertices[s][t].disp00.y=t-MatrixOfVertices[s][t].p00.y;

                MatrixOfVertices[s][t].disp01.x=s-MatrixOfVertices[s][t].p01.x;
                MatrixOfVertices[s][t].disp01.y=t-MatrixOfVertices[s][t].p01.y;

                MatrixOfVertices[s][t].disp10.x=s-MatrixOfVertices[s][t].p10.x;
                MatrixOfVertices[s][t].disp10.y=t-MatrixOfVertices[s][t].p10.y;

                MatrixOfVertices[s][t].disp11.x=s-MatrixOfVertices[s][t].p11.x;
                MatrixOfVertices[s][t].disp11.y=t-MatrixOfVertices[s][t].p11.y;

               /* e<<MatrixOfVertices[s][t].disp00.x<<'\t'<<MatrixOfVertices[s][t].disp00.y<<'\n';
                e<<MatrixOfVertices[s][t].disp01.x<<'\t'<<MatrixOfVertices[s][t].disp01.y<<'\n';
    q            e<<MatrixOfVertices[s][t].disp10.x<<'\t'<<MatrixOfVertices[s][t].disp10.y<<'\n';
                e<<MatrixOfVertices[s][t].disp11.x<<'\t'<<MatrixOfVertices[s][t].disp11.y<<'\n';
*/

                DotProduct(s,t);
                Fade(s,t);
            }
        }
    for(i=0;i<length;i=i+freq)
    {
        for(j=4;j<length;j=j+freq)
        {
            ///mai sunt puncte intre care trbeuie interpolat ca sa nu arate urat
            int r=NonLinearInterpolation(0.5);
            MatrixOfVertices[i][j].finalHeight=lerp(MatrixOfVertices[i][j-1].finalHeight,MatrixOfVertices[i][j+1].finalHeight,0.5);
        }
    }



    /*for(i=0;i<length;i++)
        for(j=0;j<length;j++)
        {
            MatrixOfVertices[i][j].finalHeight=mapFunction(MatrixOfVertices[i][j].finalHeight,min1,max1,-1,1);

        }


   /* cout<<MatrixOfVertices[2][2].q00<<endl;
        cout<<MatrixOfVertices[2][2].q01<<endl;
            cout<<MatrixOfVertices[2][2].q10<<endl;
                cout<<MatrixOfVertices[2][2].q11<<endl;
       /*for(i=0;i<128;i++)
        {
            for(j=0;j<128;j++)
            {
                g<<MatrixOfVertices[i][j].disp00.x<<'\t'<<MatrixOfVertices[i][j].disp00.y<<'\n';
                g<<MatrixOfVertices[i][j].disp01.x<<'\t'<<MatrixOfVertices[i][j].disp01.y<<'\n';
                g<<MatrixOfVertices[i][j].disp10.x<<'\t'<<MatrixOfVertices[i][j].disp10.y<<'\n';
                g<<MatrixOfVertices[i][j].disp11.x<<'\t'<<MatrixOfVertices[i][j].disp11.y<<'\n';
            }

            g<<endl;

        }*/

}

void ShowHeightMap1()
{
    int i,j;
    for(i=0;i<length;i++)
    {
        for(j=0;j<length;j++)
        {
            l<<MatrixOfVertices[i][j].finalHeight<<"f"<<' '<<","<<' ';


        }
        l<<'\n';
    }
}

void ShowHeightMap2()
{
    int i,j;
    for(i=0;i<length;i++)
    {
        for(j=0;j<length;j++)
        {
            g<<MatrixOfVertices[i][j].finalHeight<<"f"<<' '<<","<<' ';


        }
        g<<'\n';
    }
}
void ShowHeightMap3()
{
    int i,j;
    for(i=0;i<length;i++)
    {
        for(j=0;j<length;j++)
        {
            e<<MatrixOfVertices[i][j].finalHeight<<"f"<<' '<<","<<' ';


        }
        e<<'\n';
    }
}
int scalesCount;

void ResetValues()
{
    int i,j;

    for(i=0;i<129;i++)
    for(j=0;j<129;j++)
    {
        MatrixOfVertices[i][j].p00.x=0;
        MatrixOfVertices[i][j].p00.y=0;

        MatrixOfVertices[i][j].p01.x=0;
        MatrixOfVertices[i][j].p01.y=0;

        MatrixOfVertices[i][j].p10.x=0;
        MatrixOfVertices[i][j].p10.y=0;

        MatrixOfVertices[i][j].p11.x=0;
        MatrixOfVertices[i][j].p11.y=0;

        MatrixOfVertices[i][j].finalHeight=0;

        MatrixOfVertices[i][j].disp00.x=0;
        MatrixOfVertices[i][j].disp00.y=0;

        MatrixOfVertices[i][j].disp01.x=0;
        MatrixOfVertices[i][j].disp01.y=0;

        MatrixOfVertices[i][j].disp10.x=0;
        MatrixOfVertices[i][j].disp10.y=0;

        MatrixOfVertices[i][j].disp11.x=0;
        MatrixOfVertices[i][j].disp11.y=0;

        MatrixOfVertices[i][j].isLattice=false;

        MatrixOfVertices[i][j].q00=0;
        MatrixOfVertices[i][j].q01=0;
        MatrixOfVertices[i][j].q10=0;
        MatrixOfVertices[i][j].q11=0;
    }

}

int main()
{
    time_t timer=10;
    srand(time(&timer));


    f>>length;
    f>>frequency;

    AssignHashTable();
    AssignPseudGradVec(); ///assigns to all point without checking whether isLattice or not
    AssignLatticePoints(frequency);
    PickSamplePoint(frequency);
    ShowHeightMap1();
    ResetValues();


   /* f>>length;
    f>>frequency;
    AssignLatticePoints(frequency);
    PickSamplePoint(frequency);
    ShowHeightMap2();
    ResetValues();

    f>>length;
    f>>frequency;
    AssignLatticePoints(frequency);
    PickSamplePoint(frequency);
    ShowHeightMap3();
*/
///Success!
///Es ist mir gelungen! :)



    f.close();
    g.close();
    h.close();
    e.close();
    k.close();
    l.close();

    return 0;
}
