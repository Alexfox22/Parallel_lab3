#define MSMPI_NO_DEPRECATE_20
#include <iostream>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cstdlib>
#include <vector>
#include <string>
#define MPI_TAG 0
#define max_value 10


using namespace std;

struct node
{
	vector<int> neighbour;			//вектор номеров соседей
	vector<int> waight;			//вектор весов соответствующих рёбер
	int top;			//номер вершины
	int color;			//множество,которому принадлежит вершина 
	node()
	{
		color = 0;
	}
};
struct rib {
	int beg;
	int end;
	int wai;
	rib()
	{
		beg = 0;
		end = 0;
		wai = 0;
	}
};

MPI_Datatype mpi_rib;     //  Тип данных
node * graph;			    //массив вершин,в которых храним граф и множество
rib * Kr;    //исходный массив рёбер, с которым работает последовательный алгоритм (будет отсортирован...испорчен короче)
rib * KrParallel;    //исходный массив рёбер, с которым работает параллельный алгоритм (разобьет на куски, разошлет и тд)
rib * resPosl;    //результат
rib * resParallel;    //результат
int mpi_rank, mpi_size;            //переменные кол-ва процессов и ранга процесса
int tops;
int edges;
int part;
int RIB;




void let_me_see(double *mas, int len) 
{   
	
	for (int i = 0; i < len; i++)
	{
		cout << mas[i] << "  ";
	}
	cout << endl;
}
void outp(int **mas, int tops)
{
	for (int j = 0;j < tops;j++)
	{
		for (int i = 0;i < tops;i++)
		{
			if (mas[i][j] > INT_MAX - 10)
				cout << "0  ";
			else {
				cout << mas[j][i];
				if ((mas[i][j] - 9) > 0)
					cout << " ";
				else cout << "  ";
			}
		}
		cout << endl;
	}
	cout << endl;
}
void outprib(int RIB, rib *Kr)
{
	for (int o = 0;o<RIB;o++)
	{
		cout << "RIB: " << Kr[o].beg << "---" << Kr[o].end << " | " << Kr[o].wai << endl;
	}
}
bool compareres(rib *first, rib *second,int RIB)
{
	bool res = true;
	for (int i = 0;i < RIB;i++)
	{
		if ((first[i].beg != second[i].beg) || (first[i].end != second[i].end)
			|| (first[i].wai != second[i].wai))
		{
			res = false;
			break;
		}
	}
	return res;

}
void make_null(int** matrix, int len)
{
	for (int i = 0;i< len;i++)
		for (int j = 0;j < len;j++)
			matrix[i][j] = 0;
}
void generatematrix(int** mat, int len, int edges, int first_value, int last_value)
{
	bool ind = true;
	int i = 0;
	//cout << edges << "---EDGES" << endl;
	while (i<len - 2)
	{
		int x = i + 1 + rand() % (len - i - 1);
		int value = first_value + rand() % (last_value - first_value + 1);
		if ((x != i) && (mat[i][x] == 0))
		{

			mat[i][x] = value;
			mat[x][i] = value;
			i++;
		}
	}

	while (i<edges)
	{
		int x = rand() % (len);
		int y = rand() % (len);
		int value = first_value + rand() % (last_value - first_value + 1);
		if ((x != y) && (mat[x][y] == 0))
		{
			mat[x][y] = value;
			mat[y][x] = value;
			//cout << "dobavil 1 k " << i << endl;
			i++;
		}
	}

		//cout<<endl << i << endl;
	for (int i = 0;i < len;i++)
	{
		int countnn = 0;
		for (int j = 0;j < len;j++)
		{
			if (mat[i][j] > 0)
			{
				countnn++;
			}
		}
		//cout << countnn << "---COUNTNN" << endl;
		if ((countnn == 0) && (i != len - 1))
		{
			ind = false;
	//		cout << "false " << i << endl;
		}
	}
}
void generategraph(node* gr, int** mat, int tops, int edges)
{
	for (int i = 0;i < tops;i++)
	{
		gr[i].top = i + 1;
		for (int j = 0;j < tops;j++)
		{
			if (mat[i][j]>0)
			{
				gr[i].neighbour.push_back(j + 1);
				gr[i].waight.push_back(mat[i][j]);
			}
		}
	}
}
void make_matrix_wei(node *gr, int tops, int **mas)
{
	for (int i = 0;i < tops;i++)
	{
		for (int j = 0;j < tops;j++)
		{
			mas[i][j] = INT_MAX - 1;
		}
	}
	for (int i = 0;i<tops;i++)
	{
		for (int j = 0;j < gr[i].neighbour.size();j++)
		{
			mas[i][gr[i].neighbour[j] - 1] = gr[i].waight[j];
		}
	}
}
void getgraph(node *gr, int tops)
{

	for (int i = 0;i < tops;i++)
	{

		cout << "top " << gr[i].top << "  ";
		for (int j = 0;j < gr[i].neighbour.size();j++)
			cout << gr[i].neighbour[j] << ' ';
		cout << endl;

	}

}
void writegraph(node *resPrim, int tops)
{
	int checktops = 0;
	//	ofstream exfile;
	//exfile.open("outputPrim.txt");
	//string result;
	{	for (int i = 0;i < tops;i++)		//вывод
	{
		int f = 0;
		if (resPrim[i].top >= 0)
		{
			//result = result + to_string(resPrim[i].top) + "\n";
			cout << resPrim[i].top << endl;
			while (f < resPrim[i].neighbour.size())
			{
				checktops++;
				//result = result + to_string(resPrim[i].neighbour[f]) + "  ";
				cout << resPrim[i].neighbour[f] << "  ";
				f++;
			}
			//result = result + "\n";
			cout << endl;
			int fu = 0;
			while (fu < resPrim[i].waight.size())
			{
				//result = result + to_string(resPrim[i].waight[fu]) + "  ";
				cout << resPrim[i].waight[fu] << "  ";
				fu++;
			}
			//result = result + "\n";
			cout << endl << endl;
			//result = result + "\n";
		}
	}
	}
	//result = result + to_string(time / 10000.000);
	//exfile << result;
	//exfile.close();
	cout << "CHECKED----------->" << checktops << endl;
}
void fullmasofribs(rib *Kr,  int tops,node *graph)
{
	int u = 0;

	for (int j = 0;j < tops;j++)
	{
		for (int k = 0;k < graph[j].neighbour.size();k++)
		{
			/*if ((j+1 < graph[j].neighbour[k])||(j==0))
			{*/
				Kr[u].beg = j + 1;
				//cout << "***" << Kr[u].beg << "***" << endl;
				Kr[u].end = graph[j].neighbour[k];
				Kr[u].wai = graph[j].waight[k];
			//	cout << u << " RIB: " << Kr[u].beg << "---" << Kr[u].end << " | " << Kr[u].wai << endl;
				//cout << "RIBRIBRIBRIB----------" << RIB << "----------RIBRIBRIBRIB" << endl;
				u++;
			//}
		}
	}

	//cout << "----------BARRIER----------" << endl;


}
rib* algKrus(rib *Kr,int RIB,node *graph,int tops)
{
	
//	outprib(RIB, Kr);
	int i = 0;
	int countmn = 2;
	//string result1;
	rib *res = new rib[tops - 1];
	int count=0;
	graph[Kr[0].beg - 1].color = 1;
	//cout << Kr[0].beg << "__its color was changed to 1" << endl;
	for (int j = 0;j < RIB;j++)
	{
		//for(int g=0;g<RIB;g++)
		if (((graph[Kr[j].beg - 1].color == 0) && (graph[Kr[j].end - 1].color != 0))|| ((graph[Kr[j].beg - 1].color != 0) && (graph[Kr[j].end - 1].color == 0)))//||((masKr[Kr[j].beg - 1].color != masKr[Kr[j].end - 1].color)))
		{


			//result1 = result1 + to_string(Kr[j].beg) + " ---- " + to_string(Kr[j].end) + " | " + to_string(Kr[j].wai) + "\n";
		//	cout << Kr[j].beg << " ---- " << Kr[j].end << " | " << Kr[j].wai << endl;
			res[count].beg = Kr[j].beg;
			res[count].end = Kr[j].end;
			res[count].wai = Kr[j].wai;
			count++;
			if (graph[Kr[j].beg - 1].color == 0)
			{
				graph[Kr[j].beg - 1].color = graph[Kr[j].end - 1].color;

			}
			if (graph[Kr[j].end - 1].color == 0)
			{
				graph[Kr[j].end - 1].color = graph[Kr[j].beg - 1].color;

			}
		}
		if ((graph[Kr[j].beg - 1].color == 0) && (graph[Kr[j].end - 1].color == 0))
		{


		//	result1 = result1 + to_string(Kr[j].beg) + " ---- " + to_string(Kr[j].end) + " | " + to_string(Kr[j].wai) + "\n";
			//	cout << Kr[j].beg << " ---- " << Kr[j].end << " | " << Kr[j].wai << endl;
				res[count].beg = Kr[j].beg;
				res[count].end = Kr[j].end;
				res[count].wai = Kr[j].wai;
			count++;
			graph[Kr[j].beg - 1].color = countmn;

			graph[Kr[j].end - 1].color = countmn;

			countmn++;
		}
		if (((graph[Kr[j].beg - 1].color != graph[Kr[j].end - 1].color)) && ((graph[Kr[j].beg - 1].color > 0) && (graph[Kr[j].end - 1].color > 0)))
		{

			//	cout << Kr[j].beg << " ---- " << Kr[j].end << " | " << Kr[j].wai << endl;
				res[count].beg = Kr[j].beg;
				res[count].end = Kr[j].end;
				res[count].wai = Kr[j].wai;
			count++;
			//result1 = result1 + to_string(Kr[j].beg) + " ---- " + to_string(Kr[j].end) + " | " + to_string(Kr[j].wai) + "\n";
			if (graph[Kr[j].end - 1].color > graph[Kr[j].beg - 1].color)
			{
				int memory;
				//	cout << "hey, " << graph[Kr[j].end - 1].color << " > " << graph[Kr[j].beg - 1].color << endl;
				memory = graph[Kr[j].end - 1].color;
				graph[Kr[j].end - 1].color = graph[Kr[j].beg - 1].color;
				//cout << "top " << Kr[j].end << " was changed, color: " << graph[Kr[j].end - 1].color;
				//cout << endl;
				for (int y = 0;y < tops;y++)
				{
					if (graph[y].color == memory)
					{
						//		cout << "--color--" << graph[y].color;
						graph[y].color = graph[Kr[j].beg - 1].color;
						//		cout << " from " << graph[y].top << " <--- was changed " << endl;
					}
				}
			}
			else //if (graph[Kr[j].end - 1].color < graph[Kr[j].beg - 1].color)
			{
				int memory;
				//cout << "hey, " << graph[Kr[j].beg - 1].color << " > " << graph[Kr[j].end - 1].color << endl;
				memory = graph[Kr[j].beg - 1].color;
				graph[Kr[j].beg - 1].color = graph[Kr[j].end - 1].color;
				//	cout << "top " << Kr[j].beg << " was changed, color: " << graph[Kr[j].beg - 1].color;
				//	cout << endl;
				for (int y = 0;y < tops;y++)
				{
					if (graph[y].color == memory)
					{
						//cout << "--color--" << graph[y].color;
						graph[y].color = graph[Kr[j].end - 1].color;
						//	cout << " from " << graph[y].top << " <--- was changed " << endl;
					}
				}
			}
		}

		i++;
	}
	//cout<<count << endl;
	return res;
	
}
int merge(rib *ina, int lena, rib *inb, int lenb, rib *out) {
	int i, j;
	int outcount = 0;

	for (i = 0, j = 0; i<lena; i++) {
		while ((inb[j].wai < ina[i].wai) && j < lenb) {
			out[outcount++] = inb[j++];
		}
		out[outcount++] = ina[i];
	}
	while (j<lenb)
		out[outcount++] = inb[j++];

	return 0;
}
int domerge_sort(rib *a, int start, int end, rib *b) {
	if ((end - start) <= 1) return 0;

	int mid = (end + start) / 2;
	domerge_sort(a, start, mid, b);
	domerge_sort(a, mid, end, b);
	merge(&(a[start]), mid - start, &(a[mid]), end - mid, &(b[start]));
	for (int i = start; i<end; i++)
		a[i] = b[i];

	return 0;
}
int merge_sort(int n, rib *a) {
	rib *b = new rib[n];
	domerge_sort(a, 0, n, b);
	return 0;
}
void MPI_Pairwise_Exchange(int localn, rib *locala, int sendrank, int recvrank, MPI_Comm comm)
{

	
	int rank;
	rib *remote = new rib[localn];
	rib *all = new rib[2 * localn];
	const int mergetag = 1;
	const int sortedtag = 2;

	MPI_Comm_rank(comm, &rank);
	if (rank == sendrank) {
		MPI_Send(locala, localn, mpi_rib, recvrank, mergetag, MPI_COMM_WORLD);
		MPI_Recv(locala, localn, mpi_rib, recvrank, sortedtag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	else {
		MPI_Recv(remote, localn, mpi_rib, sendrank, mergetag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		merge(locala, localn, remote, localn, all);

		int theirstart = 0, mystart = localn;
		if (sendrank > rank) {
			theirstart = localn;
			mystart = 0;
		}
		MPI_Send(&(all[theirstart]), localn, mpi_rib, sendrank, sortedtag, MPI_COMM_WORLD);
		for (int i = mystart; i<mystart + localn; i++)
			locala[i - mystart] = all[i];
	}
}
int MPI_OddEven_Sort(int n, rib *a, int root, MPI_Comm comm)
{
	int rank, size, i;
	rib *local_a;

	// get rank and size of comm
	MPI_Comm_rank(comm, &rank); //&rank = address of rank
	MPI_Comm_size(comm, &size);

	local_a = (rib *)calloc(n / size, sizeof(rib));
	//if (rank == root)
		//if (n<31) let_me_see(a, n);

	// scatter the array a to local_a
	MPI_Scatter(a, n / size, mpi_rib, local_a, n / size, mpi_rib,
		root, comm);
	// sort local_a
	merge_sort(n / size, local_a);

	//odd-even part
	for (i = 1; i <= size; i++) {

		//printstat(rank, i, "before", local_a, n / size);

		if ((i + rank) % 2 == 0) {  // means i and rank have same nature
			if (rank < size - 1) {
				MPI_Pairwise_Exchange(n / size, local_a, rank, rank + 1, comm);
			}
		}
		else if (rank > 0) {
			MPI_Pairwise_Exchange(n / size, local_a, rank - 1, rank, comm);
		}

	}

	//printstat(rank, i - 1, "after", local_a, n / size);

	// gather local_a to a
	MPI_Gather(local_a, n / size, mpi_rib, a, n / size, mpi_rib,
		root, comm);

	//if (rank == root)
	//	if (n < 31) let_me_see(a, n);//printstat(rank, i, " all done ", a, n);

	return MPI_SUCCESS;
}
void compare_exchange(rib a, rib b)
{
//	cout << "here1" << endl;
	//cout << a.wai << " " << b.wai << endl;
	if (b.wai > a.wai)
	{
	//	cout << "here" << endl;
		rib tmp;
		tmp.beg = a.beg;
		tmp.end = a.end;
		tmp.wai = a.wai;
		a.beg = b.beg;
		a.end = b.end;
		a.wai = b.wai;
		b.beg = tmp.beg;
		b.end = tmp.end;
		b.wai = tmp.wai;
	}
}
void OddEvenSort(rib *A, int n) {
	//cout << "NNNNNN " << n << endl;
	for (int i = 1; i < n; i++) {
		if (i % 2 == 1) {    // нечетная итерация 
			for (int j = 0; j < n / 2 - 2; j++)
				compare_exchange(A[2 * j + 1], A[2 * j + 2]);
			if (n % 2 == 1) // сравнение последней пары при нечетном n
				compare_exchange(A[n - 2], A[n - 1]);
		}
		else    // четная итерация 
			for (int j = 1; j < n / 2 - 1; j++)
				compare_exchange(A[2 * j], A[2 * j + 1]);
	}
}

int main(int argc,char** argv)
{
	MPI_Init(&argc, &argv);            //инициализация

	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	MPI_Type_contiguous(3, MPI_INT, &mpi_rib); //создаем новый типа данных состоящий из 3 интов назовем его mpi_edge
	MPI_Type_commit(&mpi_rib); //согласовываем
	MPI_Status Status;                 //переменная статуса выполнения операции приема данных
	int tops = atoi(argv[1]);  //количество вершин

		edges = (tops*(tops - 1)) / 2;
		part = tops / mpi_size;
		double startTime, endTime;
		int first_valueRAN = 1;
		int last_valueRAN = 100;
		graph = new node[tops];
		int **matrixWEI = new int*[tops];

		for (int i = 0;i < tops;i++)
		{
			matrixWEI[i] = new int[tops];
		}
		if (mpi_rank == 0)
		{
			cout << "WELCOME" << endl;

			srand(time(0));
			make_null(matrixWEI, tops);
			generatematrix(matrixWEI, tops, edges, first_valueRAN, last_valueRAN);
		//	if (tops < 31) outp(matrixWEI, tops);
			generategraph(graph, matrixWEI, tops, edges);

			make_matrix_wei(graph, tops, matrixWEI);
			if (tops < 31) outp(matrixWEI, tops);
			//	if (tops < 31) getgraph(graph, tops);
			//	if (tops < 31) writegraph(graph, tops);
			int RIB = edges *2;
		//	if (tops < 31) cout << "RIBRIBRIBRIB----------" << RIB << "----------RIBRIBRIBRIB" << endl;
			Kr = new rib[RIB];
			KrParallel = new rib[RIB];
			fullmasofribs(Kr, tops, graph);
			fullmasofribs(KrParallel, tops, graph);
			if (tops < 6) outprib(RIB, Kr);

			//	if (tops < 31)
				//	outprib(RIB, Kr);
			resPosl = new rib[tops - 1];
			startTime = MPI_Wtime();
			//	oddEvenSorting(Kr, RIB);

			OddEvenSort(Kr, RIB);
			cout << endl;
			if (tops < 6) outprib(RIB, Kr);
			resPosl = algKrus(Kr, RIB, graph, tops);
			cout << RIB << endl;
			endTime = MPI_Wtime();
			cout << endl;
			if (tops < 31) outprib(tops - 1, resPosl);

			double t1 = endTime - startTime;
			cout << "Spended time: " << t1 << endl;

		}
		if (mpi_rank == 0)                 //если в root
		{
			for (int i = 0;i < tops;i++)
				graph[i].color = 0;
			startTime = MPI_Wtime();
			cout << "mpi_size=" << mpi_size << endl; //вывод кол-ва действующих процессов в программе
			cout << endl;
			startTime = MPI_Wtime();
		}

		MPI_OddEven_Sort(RIB, KrParallel, 0, MPI_COMM_WORLD);
		if (mpi_rank == 0)

		{
			cout << endl;
			if (tops < 31) outprib(RIB, KrParallel);
			resParallel = new rib[tops - 1];
			resParallel = algKrus(KrParallel, RIB, graph, tops);
			endTime = MPI_Wtime();
			cout << endl;
			double t2 = endTime - startTime;
			if (tops < 31) outprib(tops - 1, resPosl);
			if (compareres(resPosl, resParallel, RIB) == true)
			{
				cout << "equal" << endl;
			}
			else cout << "not equal" << endl;
			cout << "Spended time: " << t2 << endl;
			//	cout << "Time spended: " << endTime - startTime << endl;
			cout << "------------------------------------------------" << endl;
		}

		delete[] graph;
		delete[] Kr, KrParallel;
		delete[] resParallel, resPosl;
		for (int i = 0;i < tops;i++)
			delete[] matrixWEI[i];
		delete[] matrixWEI;
	
	MPI_Finalize(); //завершаем работу с MPI
	  return 0;
}
