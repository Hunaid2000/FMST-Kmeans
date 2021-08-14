#include "Dataset.h"
#include <algorithm>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <fstream>
using namespace std;

// This constructor reads 2d-datapoints from file
Dataset::Dataset(string fname, int n) {
    N = n;
    K = sqrt(n);
    clusters.resize(K);

    ifstream in;
    in.open(fname);
    double x, y;
    Point p;
    while (!in.eof()) {
        in >> x >> y;
        datapoints.push_back(Point(x, y));
    }
    in.close();
}

// This function checks all datapoints that which cluster it should belong based on Distance Formula
void Dataset::assignPoints() {
    for (int i = 0; i < N; i++) {
        double mindist = INT_MAX;
        for (int j = 0; j < K; j++) {
            double dist = datapoints[i].calDistance(clusters[j].centeriod);
            if (dist < mindist) {
                mindist = dist;
                datapoints[i].clusterID = j + 1;
            }
        }
    }
}

// This function checks Old and New centeriods
bool Dataset::compareCenters(vector<Point>& Newcenters) {
    for (int i = 0; i < K; i++) {
        if (clusters[i].centeriod != Newcenters[i])
            return false;
    }
    return true;
}

// This function applies Kmeans to Dataset and make clusters
void Dataset::Kmeans() {
    vector<Point> centers;
    while (centers.size() < K)
    {
        int rid = rand() % N;
        if (std::find(centers.begin(), centers.end(), datapoints[rid]) == centers.end())
        {
            centers.push_back(datapoints[rid]);
        }
    }

    for (int i = 0; i < centers.size(); i++) {
        clusters[i].centeriod = centers[i];
    }

    while (true)
    {
        //assigning points to clusters
        assignPoints();
        cout << "Centeriods:" << endl;
        for (int i = 0; i < centers.size(); i++) {
            cout << clusters[i].centeriod << endl;
        }
        cout << endl;
        cout << "-----Points with ClusterID-----" << endl;
        for (int i = 0; i < N; i++) {
            cout << datapoints[i] << ": ID = " << datapoints[i].clusterID << endl;
        }
        cout << endl;

        //recompute means of centeriods
        vector<double> total(K, 0);
        vector<double> sumx(K, 0);
        vector<double> sumy(K, 0);
        for (int i = 0; i < N; i++)
        {
            int cID = datapoints[i].clusterID;
            total[cID - 1]++;
            sumx[cID - 1] += datapoints[i].x;
            sumy[cID - 1] += datapoints[i].y;
        }

        //assign new centers
        vector<Point> Newcenters(K);
        for (int i = 0; i < K; i++)
        {
            Point temp;
            temp.x = sumx[i] / total[i];
            temp.y = sumy[i] / total[i];
            Newcenters[i] = temp;
        }

        // check if there is change in centeriods
        if (compareCenters(Newcenters) == 1) {
            break;
        }
        else {
            for (int i = 0; i < K; i++)
                clusters[i].centeriod = Newcenters[i];
        }

    }

    // Assigning datapoints to each cluster
    for (int i = 0; i < N; i++) {
        int ind = datapoints[i].clusterID;
        clusters[ind - 1].clusterpoints.push_back(datapoints[i]);
    }

    cout << "Centeriods:" << endl;
    for (int i = 0; i < K; i++) {
        cout << clusters[i].centeriod << endl;
    }
    cout << endl;

}

// This function makes MST of each clusters by calling function from cluster struct
void Dataset::makeMST() {
    for (int i = 0; i < K; i++) {
        clusters[i].makeClusterMST();
        cout << "MST of ClusterID : " << i + 1 << endl;
        clusters[i].printMST();
        cout << endl;
    }

}

// Divide and Conquer Using K-means 
void Dataset::DAC() {
    Kmeans();
    makeMST();
}

// This function calculates the midpoints of MST of centers
Point Dataset::calMSTCentralMidPoint(Point p1, Point p2) {
    Point p;
    p.x = (p1.x + p2.x) / 2;
    p.y = (p1.y + p2.y) / 2;
    return p;
}

// Function to find min weight vertex
int Dataset::findMinWeightVertex(vector<bool>& visited, vector<double>& key, int size) {
    int mindist = INT_MAX, index;
    for (int i = 0; i < size; i++) {
        if (key[i] < mindist && visited[i] == false) {
            mindist = key[i];
            index = i;
        }
    }
    return index;
}

// This function is Combine Algorithm function that combines the MSTs of clusters based on MST of center
void Dataset::CA(vector<ConnectingEdge>& MST)
{
    for (int i = 0; i < K; i++)
    {
        clusters[i].centeriod.clusterID = i + 1;
    }

    // Making a graph matrix for the representation of edges
    vector<vector<double>> graphMatrix;
    for (int i = 0; i < K; i++) {
        graphMatrix.push_back(vector<double>());
        for (int j = 0; j < K; j++) {
            graphMatrix[i].push_back(clusters[i].centeriod.calDistance(clusters[j].centeriod));
        }
    }

    // Apply Prims Algorithm to find MST of centers
    MSTCentral.clear();
    MSTCentral.resize(K);
    vector<double> key(K);
    vector<bool> visited(K);
    for (int i = 0; i < K; i++) {
        key[i] = INT_MAX;
        visited[i] = false;
    }
    key[0] = 0;
    MSTCentral[0] = Point(INT_MIN, INT_MIN);

    for (int i = 0; i < K - 1; i++) {
        int m = findMinWeightVertex(visited, key, K);
        visited[m] = true;
        for (int n = 0; n < K; n++) {
            if (graphMatrix[m][n] != 0 && visited[n] == false && graphMatrix[m][n] < key[n]) {
                MSTCentral[n] = clusters[m].centeriod;
                key[n] = graphMatrix[m][n];
            }
        }
    }

    printMSTCentral();
    cout << endl;

    // Finding the edge between two MSTs
    ConnectEdges.clear();
    for (int i = 1; i < K; i++)
    {
        ConnectEdges.push_back(DCE(MSTCentral[i], clusters[i].centeriod));
    }

    // Combining the MSTs of all Clusters
    MST.resize(N);
    int n = 1;
    for (int i = 0; i < K; i++)
    {
        for (int j = 1; j < clusters[i].MST.size(); j++) {
            MST[n].PointA = clusters[i].MST[j];
            MST[n].PointB = clusters[i].clusterpoints[j];
            n++;
        }
    }
    for (int i = 0; i < ConnectEdges.size(); i++) {
        MST[n].PointA = ConnectEdges[i].PointA;
        MST[n].PointB = ConnectEdges[i].PointB;
        n++;
    }

}

// This function Detects the Connecting Edge between two MSTs
ConnectingEdge Dataset::DCE(Point p1, Point p2)
{
    int min = INT_MAX;
    Point a, b;

    // finding point a in 1st Cluster that is closest to centeriod of 2nd Cluster
    for (int i = 0; i < datapoints.size(); i++) {
        if (datapoints[i].clusterID == p2.clusterID) {
            double dist = datapoints[i].calDistance(p1);
            if (dist < min) {
                min = dist;
                b = datapoints[i];
            }
        }
    }

    // finding point b in 2nd Cluster that is closest to centeriod of 1st Cluster
    min = INT_MAX;
    for (int i = 0; i < datapoints.size(); i++) {
        if (datapoints[i].clusterID == p1.clusterID) {
            double dist = datapoints[i].calDistance(p2);
            if (dist < min) {
                min = dist;
                a = datapoints[i];
            }
        }
    }

    ConnectingEdge Edge;
    Edge.PointA = a;
    Edge.PointB = b;
    return Edge;
}

// This function finds Secondary Approximate MST
void Dataset::SAM() {
    vector<Point> newcenters;
    for (int i = 1; i < K; i++) {
        newcenters.push_back(calMSTCentralMidPoint(MSTCentral[i], clusters[i].centeriod));
    }
    // Partitioning the Dataset based on Refinement Stage
    clusters.clear();
    clusters.resize(K - 1);

    cout << "Centeriods:" << endl;
    for (int i = 0; i < clusters.size(); i++) {
        clusters[i].centeriod = newcenters[i];
        cout << clusters[i].centeriod << endl;
    }
    cout << endl;

    for (int i = 0; i < datapoints.size(); i++) {
        datapoints[i].clusterID = 0;
    }

    // Assiging Points to new clusters
    this->K = K - 1;
    assignPoints();
    for (int i = 0; i < N; i++) {
        int ind = datapoints[i].clusterID;
        clusters[ind - 1].clusterpoints.push_back(datapoints[i]);
    }

    // Making MSTs of each cluster and combining them
    cout << endl;
    makeMST();
    CA(MST2);
    printMST2();
}


void Dataset::printMST1() {
    cout << "---------MST1---------" << endl;
    for (int i = 1; i < MST1.size(); i++)
    {
        cout << MST1[i].PointA << " - " << MST1[i].PointB << endl;
    }
}

void Dataset::printMST2() {
    cout << "---------MST2---------" << endl;
    for (int i = 1; i < MST2.size(); i++)
    {
        cout << MST2[i].PointA << " - " << MST2[i].PointB << endl;
    }
}

void Dataset::printFinalMST() {
    cout << "---------FinalMST---------" << endl;
    for (int i = 1; i < FinalMST.size(); i++)
    {
        cout << FinalMST[i] << " - " << datapoints[i] << endl;
    }
}

void Dataset::printMSTCentral() {
    cout << "-------MSTofCenters-------" << endl;
    for (int i = 1; i < K; i++)
        cout << MSTCentral[i].clusterID << " - " << clusters[i].centeriod.clusterID << endl;
}

// A utility function to find index of point
int Dataset::findind(Point p) {
    for (int i = 0; i < datapoints.size(); i++) {
        if (p.x == datapoints[i].x && p.y == datapoints[i].y)
            return i;
    }
    return -1;
}

// Fast MST
void Dataset::FMST() {
    cout << "\n-------Divide And Conquer Stage-------" << endl;
    DAC();
    CA(MST1);
    printMST1();
    cout << "\n-------Refinement Stage-------" << endl;
    SAM();

    // Making a graph matrix for the representation of edges
    vector<vector<double>> graphMatrix(N, vector<double>(N, 0));
    for (int j = 1; j < N; j++) {
        int m = findind(MST1[j].PointA);
        int n = findind(MST1[j].PointB);
        graphMatrix[m][n] = MST1[j].PointA.calDistance(MST1[j].PointB);
        graphMatrix[n][m] = MST1[j].PointA.calDistance(MST1[j].PointB);
    }
    for (int j = 1; j < N; j++) {
        int m = findind(MST2[j].PointA);
        int n = findind(MST2[j].PointB);
        graphMatrix[m][n] = MST2[j].PointA.calDistance(MST2[j].PointB);
        graphMatrix[n][m] = MST2[j].PointA.calDistance(MST2[j].PointB);
    }

    // Apply Prims Algorithm to final Approximate MST
    FinalMST.resize(N);
    vector<double> key(N);
    vector<bool> visited(N);
    for (int i = 0; i < N; i++) {
        key[i] = INT_MAX;
        visited[i] = false;
    }
    key[0] = 0;
    FinalMST[0] = Point(INT_MIN, INT_MIN);

    for (int i = 0; i < N - 1; i++) {
        int m = findMinWeightVertex(visited, key, N);
        visited[m] = true;
        for (int n = 0; n < N; n++) {
            if (graphMatrix[m][n] != 0 && visited[n] == false && graphMatrix[m][n] < key[n]) {
                FinalMST[n] = datapoints[m];
                key[n] = graphMatrix[m][n];
            }
        }
    }
    cout << endl;
    printFinalMST();

}
