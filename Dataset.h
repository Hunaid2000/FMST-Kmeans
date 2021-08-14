#pragma once
#include <iostream>
#include <vector>
using namespace std;

struct Point {
    double x, y;
    int clusterID;

    Point() {
        x = 0;
        y = 0;
        clusterID = 0;
    }

    Point(double x, double y) {
        this->x = x;
        this->y = y;
        clusterID = 0;
    }

    double calDistance(Point p) {
        double dist = pow(abs(p.x - this->x), 2) + pow(abs(p.y - this->y), 2);
        return sqrt(dist);
    }

    Point& operator= (const Point& p) {
        x = p.x;
        y = p.y;
        clusterID = p.clusterID;
        return *this;
    }

    bool operator== (const Point& p)
    {
        return (this->x == p.x && this->y == p.y && this->clusterID == p.clusterID);
    }

    bool operator!= (const Point& p)
    {
        return !(*this == p);
    }

    friend ostream& operator<<(ostream& os, const Point& pt) {
        os << "(" << pt.x << "," << pt.y << ")";
        return os;
    }
};


struct Cluster {
    vector <Point> MST;
    Point centeriod;
    vector <Point> clusterpoints;

    // This function returns index of minimum weight vertex edge
    int findMinWeightVertex(vector<bool>& visited, vector<double>& key) {
        int mindist = INT_MAX, index;
        for (int i = 0; i < clusterpoints.size(); i++) {
            if (key[i] < mindist && visited[i] == false) {
                mindist = key[i];
                index = i;
            }
        }
        return index;
    }

    // This function creates MST of a cluster bases on Prims Algorithm
    void makeClusterMST() {
        // Making a graph matrix for the representation of edges
        vector<vector<double>> graphMatrix;
        for (int i = 0; i < clusterpoints.size(); i++) {
            graphMatrix.push_back(vector<double>());
            for (int j = 0; j < clusterpoints.size(); j++) {
                graphMatrix[i].push_back(clusterpoints[i].calDistance(clusterpoints[j]));
            }
        }

        // Apply Prims Algorithm to find MST of Cluster
        MST.resize(clusterpoints.size());
        vector<double> key(clusterpoints.size());
        vector<bool> visited(clusterpoints.size());
        for (int i = 0; i < clusterpoints.size(); i++) {
            key[i] = INT_MAX;
            visited[i] = false;
        }
        key[0] = 0;
        MST[0] = Point(INT_MIN, INT_MIN);

        for (int i = 0; i < clusterpoints.size() - 1; i++) {
            int m = findMinWeightVertex(visited, key);
            visited[m] = true;
            for (int n = 0; n < clusterpoints.size(); n++) {
                if (graphMatrix[m][n] != 0 && visited[n] == false && graphMatrix[m][n] < key[n]) {
                    MST[n] = clusterpoints[m];
                    key[n] = graphMatrix[m][n];
                }
            }
        }

    }

    void printMST() {
        for (int i = 1; i < clusterpoints.size(); i++)
            cout << MST[i] << " - " << clusterpoints[i] << endl;
    }

};


struct ConnectingEdge
{
    Point PointA;
    Point PointB;
};


class Dataset {
    int N, K;
    vector <Point> datapoints;
    vector <Cluster> clusters;
    vector <Point> MSTCentral;
    vector <ConnectingEdge> ConnectEdges;
    vector <ConnectingEdge> MST1;
    vector <ConnectingEdge> MST2;
    vector <Point> FinalMST;
public:
    // This constructor reads 2d-datapoints from file
    Dataset(string fname, int n);

    // This function checks all datapoints that which cluster it should belong based on Distance Formula
    void assignPoints();

    // This function checks Old and New centeriods
    bool compareCenters(vector<Point>& Newcenters);

    // This function applies Kmeans to Dataset and make clusters
    void Kmeans();

    // This function makes MST of each clusters by calling function from cluster struct
    void makeMST();

    // Divide and Conquer Using K-means 
    void DAC();

    // This function calculates the midpoints of MST of centers
    Point calMSTCentralMidPoint(Point p1, Point p2);

    // Function to find min weight vertex
    int findMinWeightVertex(vector<bool>& visited, vector<double>& key, int size);

    // This function is Combine Algorithm function that combines the MSTs of clusters based on MST of center
    void CA(vector<ConnectingEdge>& MST);

    // This function Detects the Connecting Edge between two MSTs
    ConnectingEdge DCE(Point p1, Point p2);

    // This function finds Secondary Approximate MST
    void SAM();

    void printMST1();

    void printMST2();

    void printFinalMST();

    void printMSTCentral();

    // A utility function to find index of point
    int findind(Point p);

    // Fast MST
    void FMST();

};