#pragma once
#include "Point.h"
#include "QubeKe.h"

ofstream logFile("log.txt");
const double PI = 3.14159265;
double *xNet, *yNet, *zNet;
int nX, nY, nZ;
QubeKe*qubes;
int nQubes;
vector<double> pTheoriticAvg;
vector<double> gravity_pr; 
vector<double> gravitySolution;
double **A, *x, *b;
double alpha;
int *Jm;
vector<Point> receivers;

void readAxis(ifstream &fileNet, set<double> &mas) {
	double firstElement, startPosition, endPosition, position, lastPosition, stepReal;;
	int i, numberOfInterval;
	fileNet >> firstElement >> numberOfInterval;
	mas.insert(firstElement);
	if (numberOfInterval < 1) {
		return;
	}
	double*intervals = new double[numberOfInterval + 1];
	double*sizeOfSteps = new double[numberOfInterval];
	double*mnojiteli = new double[numberOfInterval];
	int*napravlenie = new int[numberOfInterval];
	intervals[0] = firstElement;
	for (i = 1; i < numberOfInterval + 1; i++)
	{
		fileNet >> intervals[i];
		mas.insert(intervals[i]);
	}
	for (i = 0; i < numberOfInterval; i++)
	{
		fileNet >> sizeOfSteps[i];
	}
	for (i = 0; i < numberOfInterval; i++)
	{
		fileNet >> mnojiteli[i];
		if (mnojiteli[i] < 1)
		{
			cout << "Ошибка описания сетки по оси. Коэфициент разрядки должен быть больше или равен 1";
			system("pause");
			exit(1);
		}
	}
	for (i = 0; i < numberOfInterval; i++)
	{
		fileNet >> napravlenie[i];
	}
	for (i = 0; i < numberOfInterval; i++)
	{
		if (napravlenie[i] == -1)
		{
			startPosition = intervals[i + 1];
			endPosition = intervals[i];
		}
		else
		{
			startPosition = intervals[i];
			endPosition = intervals[i + 1];
		}
		stepReal = 0;
		lastPosition = startPosition;
		position = startPosition + napravlenie[i] * sizeOfSteps[i];
		while (position*napravlenie[i] < endPosition*napravlenie[i])
		{
			mas.insert(position);
			stepReal = fabs(lastPosition - position)*mnojiteli[i];
			lastPosition = position;
			position += napravlenie[i] * stepReal;
		}
		//if (fabs(lastPosition - endPosition) < sizeOfSteps[i] && lastPosition != startPosition)
		//{
		//	mas.erase(mas.find(lastPosition));
		//}
	}
}

double* setToArray(set<double> setV, int &n) {
	n = setV.size();
	double *array = new double[n];
	int i = 0;
	for (set<double>::const_iterator it = setV.begin(); it != setV.end(); it++, i++)
	{
		array[i] = *it;
		logFile << array[i] << " ";
	}
	return array;
}

void input(string taskFile, string gridFile) {
	ifstream grid(gridFile);
	set<double> x, y, z;
	readAxis(grid, x);
	readAxis(grid, y);
	readAxis(grid, z);

	logFile << endl << "X:" << endl;
	xNet = setToArray(x, nX);
	logFile << endl << "Y:" << endl;
	yNet = setToArray(y, nY);
	logFile << endl << "Z:" << endl;
	zNet = setToArray(z, nZ);
	logFile << endl;

	ifstream taskDeskriptor(taskFile);
	taskDeskriptor >> alpha;
	logFile << "alpha = " << alpha << endl;
}

void createReceivers(string receiversFile) {
	ifstream receiversStream(receiversFile);
	set<double> x, y,z;
	readAxis(receiversStream, x);
	readAxis(receiversStream, y);
	readAxis(receiversStream, z);

	logFile << "Receivers:" << endl;
	for each (double zValue in z)
	{
		for each (double yValue in y)
		{
			for each (double xValue in x)
			{
				Point receiver = Point(xValue, yValue, zValue);
				receivers.push_back(receiver);
				logFile << "\t" << receiver << endl;
			}
		}
	}

	logFile << endl;
}

double pTheoritic(int iQube) {
	//if (iQube == 5 || iQube == 6 || iQube == 9 || iQube == 10)
	//	return 1;
	//return 0;

	//if (iQube == 23 || iQube == 24 || iQube == 32 || iQube == 33)
	//	return 1;
	//return 0;

	return 1;
}

void createQubes() {
	nQubes = (nX - 1)*(nY - 1)*(nZ - 1);
	qubes = new QubeKe[nQubes];
	pTheoriticAvg.resize(nQubes);
	gravity_pr.resize(receivers.size());
	int iQube = 0;
	logFile << "Qubes:" << endl;
	for (size_t k = 0; k < nZ - 1; k++)
	{
		for (size_t j = 0; j < nY - 1; j++)
		{
			for (size_t i = 0; i < nX - 1; i++)
			{
				qubes[iQube].initQube(xNet[i], yNet[j], zNet[k], xNet[i + 1], yNet[j + 1], zNet[k + 1]);
				logFile << qubes[iQube] << endl;
				pTheoriticAvg[iQube] = pTheoritic(iQube);
				iQube++;
			}
		}
	}
}

double length(Point p1, Point p2) {
	return sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y) + (p1.z - p2.z)*(p1.z - p2.z));
}

double calcByGauss(int iReceiver, int iQube) {
	double hx = qubes[iQube].nodes[1].x - qubes[iQube].nodes[0].x;
	double hy = qubes[iQube].nodes[2].y - qubes[iQube].nodes[0].y;
	double hz = qubes[iQube].nodes[4].z - qubes[iQube].nodes[0].z;
	Point startPoint = qubes[iQube].nodes[0];
	double tKoef = sqrt(3. / 5.);
	double integrationPoints[3] = {
		0.5,
		(1 + tKoef) / 2,
		(1 - tKoef) / 2
	};
	double tauKoefs[3] = {
		8. / 9.,
		5. / 9.,
		5. / 9.
	};
	double sum = 0;
	for (size_t i = 0; i < 3; i++)
	{
		double xCoord = integrationPoints[i] * hx;
		for (size_t j = 0; j < 3; j++)
		{
			double yCoord = integrationPoints[j] * hy;
			for (size_t k = 0; k < 3; k++)
			{
				double zCoord = integrationPoints[k] * hz;

				double mes = qubes[iQube].calcVolume();
				//Point M1 = qubes[iQube].calcCenter();
				Point M1(startPoint.x + xCoord, startPoint.y + yCoord, startPoint.z + zCoord);
				Point M2 = receivers[iReceiver];

				Point dp = M1 - M2;
				double r3 = pow(length(M1, M2), 3);
				double koef = mes / (4.0 * PI * r3);

				sum += tauKoefs[i] * tauKoefs[j] * tauKoefs[k] * koef * dp.z;
			}
		}
	}
	return sum;
}

double calcDeltaGk_z(int iReceiver, int iQube)
{
	double gaussVal =  calcByGauss(iReceiver, iQube);

	double mes = qubes[iQube].calcVolume();
	Point M1 = qubes[iQube].calcCenter();
	Point M2 = receivers[iReceiver];

	Point dp = M1 - M2;
	double r3 = pow(length(M1, M2), 3);
	double koef = mes / (4.0 * PI * r3);
	double result = (koef * dp.z);

	return result;
}

void calcGz(vector<double> &gravity, double *power) {
	for (size_t i = 0; i < receivers.size(); i++)
	{
		double sum = 0;
		for (size_t iQube = 0; iQube < nQubes; iQube++)
		{
			double val = calcDeltaGk_z(i, iQube);
			sum += power[iQube] * val;
		}
		gravity[i] = sum;
	}
}

void buildMatrix(int n) {
	A = new double*[n];
	b = new double[n];
	x = new double[n];
	Jm = new int[n];

	for (int i = 0; i < n; i++) {
		Jm[i] = i;
	}

	for (size_t i = 0; i < n; i++)
	{
		b[i] = x[i] = 0;
		A[i] = new double[n];
		for (size_t j = 0; j < n; j++)
		{
			A[i][j] = 0;
		}
	}

	for (int i = 0; i < receivers.size(); i++)
	{
		for (int q = 0; q < n; q++)
		{
			double dgz_iq = calcDeltaGk_z(i, q);
			for (int s = 0; s < q; s++)
			{
				double dgz_is = calcDeltaGk_z(i, s);
				double val = dgz_is * dgz_iq;

				A[q][s] += val;
				A[s][q] += val;
			}
			A[q][q] += (dgz_iq * dgz_iq);
			//logFile << setw(6) << i << setw(6) << q << setw(18) << dgz_iq << endl;
			b[q] += (dgz_iq * gravity_pr[i]);
		}
	}

	for (int q = 0; q < n; q++)
		A[q][q] += alpha;

	logFile << "A:" << endl;
	for (size_t i = 0; i < n; i++)
	{
		for (size_t j = 0; j < n; j++)
		{
			logFile << setw(10) << A[i][j] << " ";
		}
		logFile << endl;
	}
	logFile << endl << "b:" << endl;
	for (size_t i = 0; i < n; i++)
	{
		logFile << b[i] << " ";
	}
}

void search_main_elem(int i)
{
	int max = i;
	int n = nQubes;
	for (int k = i + 1; k<n; k++)
		if (fabs(A[i][Jm[k]]) > fabs(A[i][Jm[max]])) max = k;

	if (max != i)
	{
		int buf = Jm[i];
		Jm[i] = Jm[max];
		Jm[max] = buf;
	}
}

void solveSLAU(int n)
{
	for (int i = 0; i<n; i++)
	{
		x[i] = 0;

		search_main_elem(i);
		double diag = A[i][Jm[i]];

		for (int k = i + 1; k<n; k++)
		{
			double m = A[k][Jm[i]] / diag;
			A[k][Jm[i]] = 0;
			for (int j = i + 1; j<n; j++)
				A[k][Jm[j]] -= (A[i][Jm[j]] * m);
			b[k] -= (b[i] * m);
		}

	}

	for (int i = n - 1; i >= 0; i--)
	{
		double sum = 0;
		for (int j = i + 1; j<n; j++)
			sum += A[i][Jm[j]] * x[Jm[j]];

		x[Jm[i]] = (b[i] - sum) / A[i][Jm[i]];
	}
}

void output() {
	ofstream outP("outP.txt");
	outP << "Qube" << setw(18) << "p_th" << setw(18) << "p_solution" << endl;

	for (int i = 0; i < nQubes; i++)
	{
		outP << setw(5) << i <<setw(18)<< pTheoriticAvg[i] <<setw(18) << x[i] << endl;
	}

	ofstream outGz("outGz.txt");
	outGz << "Receiver" << setw(18) << "gz_p_th" << setw(18) << "gz_p_solution" << setw(18) << "dGz/gz_p_th" << endl;
	double sum = 0;
	for (int i = 0; i < receivers.size(); i++)
	{
		double res = fabs(gravitySolution[i] - gravity_pr[i]);
		outGz << setw(5) << i << setw(18) << gravity_pr[i] << setw(18) << gravitySolution[i] << setw(18)
			<< res / gravity_pr[i] << endl;
		sum += res*res;
	}
	outGz << "Minimization = " << sum << endl;
}

double GZinput(Point point) {
	double sigma = 1500;
	double matOj = 2500;
	double normal_disperce = 1 / (sigma*sqrt(2 * PI))*exp(-((point.x - matOj)*(point.x - matOj)) / (2 * sigma*sigma));
	return 160000*normal_disperce;
}

void GZfromFuncAndPrintInFile() {
	logFile << "GZ" << endl;
	for (int i = 0; i < receivers.size(); i++)
	{
		double val = GZinput(receivers[i]);
		logFile << setw(8) << i << setw(18) << receivers[i] << setw(18) << val << endl;
		gravity_pr[i] = val;
	}
	logFile << endl;
}

void main() {
	input("task.txt","grid.txt");
	createReceivers("receivers.txt");
	createQubes();
	calcGz(gravity_pr, pTheoriticAvg.data());

	//GZfromFuncAndPrintInFile();

	buildMatrix(nQubes);
	solveSLAU(nQubes);

	gravitySolution.resize(receivers.size());
	calcGz(gravitySolution, x);
	output();
}