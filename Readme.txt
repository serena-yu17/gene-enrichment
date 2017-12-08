This program is called in other programs to select a sublist of genes taken from a large gene network (10k nodes) and graph their interaction through related genes. 

For example, searching for USP2, TP53 will generate: TP53 -> MDM4 -> USP2 (MDM4 is added as an intermediate gene).



The C ABI is:

extern "C" __declspec(dllexport) void enrichGenes(int** vertex, double** coordinates, int* vcount, int** edge, int* ecount, const int data[], int nData, const int genelist[], int nlist);

	int data[]:input the total list of pairs of gene id in the form of integers. e.g. a database [TP53->MDM4, AKT1-> MDM2], confidence 3, 2 will be numbered in the outside code into <data[342, 4324, 3, 5646, 13, 2]>, which represents 2 relation pairs of TP53->MDM4 of confidence 3 and AKT1->MDM2 of confidence 2.
	 
	int nData: input. Number of groups of <data[]>. It should be length of <data[]> / 3. If incorrectly set, it will cause a memory read outbond.

	int genelist[]: input all genes to be queried. e.g. a query [TP53, MDM2] will be numbered in the outside code into <genelist[342, 13]>.

	int nlist: input the length of <genelist[]>. If incorrectly set, it will cause a memory read outbond.


	int** vertex: output. Initial value should be NULL. Will allocate memory for gene id numbers to output an int array as int* vertex[].

	int* vcount: output. Number of vertices. Sould be the same as the length of <vertex>.

	double** coordinates: output. Initial value should be NULL. Will allocate memory for gene id numbers to output an int array as double* coordinates[]. Its length is 2 * <vcount>, and each pair of 2 double values represent the x and y coordinates of the vertices in <vertex>.

	int* ecount: output. Number of edges, should be <leng of <edge>> /2.

	int** edge:  output. Initial value should be NULL. Will allocate memory for gene id numbers to output an int array as int* edge[]. Its length is 2 * <ecount>, and each pair of 2 int values represent a directed edge between the 2 gene ids as integers.






extern "C" __declspec(dllexport) void freeArray(void* arr);
	Release a pointer allocated by this same program (i.e. allocated by the C++ new operator). Use for <vertex>, <coordinates>, <edge> after finishing using their data.
	


