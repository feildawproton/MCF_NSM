#include <stdlib.h>
#include <stdio.h>

typedef enum { MCF_FALSE = 0, MCF_TRUE = 1 } MCFBool;

typedef enum { A, B } MCFOptions;

typedef struct ConvLinkCandidate
{
	MCFBool Converged;
	unsigned int LinkCandidate;
	float CandidateValue;
}ConvLinkCandidate;

typedef struct PhaseIReturn
{
	unsigned int artN;
	unsigned int artL;
	unsigned int* artNi_l;
	unsigned int* artNj_l;
	MCFBool* artBasis_l;
	float* artCij_l;
	float* artXij_l;
	float* artB_n;
}PhaseIReturn;

float cost(const unsigned int L, const float* Cij_l, const float* Xij_l)
{
	float sum = 0;
	for (unsigned link = 0; link < L; link++)
	{
		sum += Cij_l[link] * Xij_l[link];
	}
	return sum;
}

MCFBool compW_Terminate(unsigned iter, const unsigned int N, const MCFBool* Wcalced_n)
{
	MCFBool Terminate = MCF_TRUE;
	for (unsigned int n = 0; n < N; n++)
	{
		if (Wcalced_n[n] == MCF_FALSE)
			Terminate = MCF_FALSE;
	}
	if (iter >= N)
	{
		Terminate = MCF_TRUE;
	}
	return Terminate;
}

//need to track if link is active
float* compW(const unsigned int N, const unsigned int L, const MCFBool* Basis_l, const unsigned int* Ni_l, const unsigned int* Nj_l, const float* Cij_l)
{
	//init array of w (W). use N as the size as they go with each node
	float* W_n = malloc(N * sizeof(float));
	MCFBool* Wcalced_n = malloc(N * sizeof(float));
	for (unsigned int n = 0; n < N; n++)
	{
		W_n[n] = 0.0f;
		if (n == N - 1)
			Wcalced_n[n] = MCF_TRUE;
		else
			Wcalced_n[n] = MCF_FALSE;
	}

	unsigned int compW_Iter = 0;
	//go until we are happy and all w's are calculated
	while (compW_Terminate(compW_Iter, N, Wcalced_n) == MCF_FALSE)
	{
		MCFBool Terminate = compW_Terminate(compW_Iter, N, Wcalced_n);
		//check the links
		for (unsigned int l = 0; l < L; l++)
		{
			//only for the bases
			if (Basis_l[l] == MCF_TRUE)
			{
				//check both combinations
				if (Wcalced_n[Nj_l[l]] == MCF_TRUE && Wcalced_n[Ni_l[l]] == MCF_FALSE)
				{
					W_n[Ni_l[l]] = Cij_l[l] + W_n[Nj_l[l]];
					Wcalced_n[Ni_l[l]] = MCF_TRUE;
				}
				else if (Wcalced_n[Ni_l[l]] == MCF_TRUE && Wcalced_n[Nj_l[l]] == MCF_FALSE)
				{
					W_n[Nj_l[l]] = W_n[Ni_l[l]] - Cij_l[l];
					Wcalced_n[Nj_l[l]] = MCF_TRUE;
				}
				//otherwise give up and go home
			}
		}

		compW_Iter++;
		if (compW_Iter > N)
		{
			printf("\nSomething went wrong inside of compW. Please exit.");
			getchar;  //hack
		}
	}

	printf("\n\nPrinting w's:");
	for (unsigned int n = 0; n < N; n++)
	{
		printf("\nNode %i, %f", n, W_n[n]);
	}

	free(Wcalced_n);

	return W_n;
}

//if there is a tie pick the first one
//stack pointer?  Don't try to free this because we never alloc!
ConvLinkCandidate* create_ConvLinkCandidate(const unsigned int N, const unsigned int L, const MCFBool* Basis_l, const unsigned int* Ni_l, const unsigned int* Nj_l, const float* W_n, const float* Cij_l)
{
	//let A = L - B (total number of links minus the number of basis)
	float Max = 0;
	unsigned int M = 0;
	unsigned int NegCount = 0;
	unsigned int BasisCount = 0;
	MCFBool AlreadyCalced = MCF_FALSE;  //ADDED
	printf("\n\nPrinting zij_cij: ");
	for (unsigned int l = 0; l < L; l++)
	{
		//looking for not basis
		if (Basis_l[l] == MCF_FALSE)
		{
			BasisCount++;

			float zij_cij = W_n[Ni_l[l]] - W_n[Nj_l[l]] - Cij_l[l];
			printf("\nLink %i, ni_l %i, nj_l %i, cij_l %f, and zij_cij value %f", l, Ni_l[l], Nj_l[l], Cij_l[l], zij_cij);

			if (zij_cij > Max)
			{
				Max = zij_cij;
				M = l;
				AlreadyCalced = MCF_TRUE;
			}
			////ADDED THIS MIGHT NEED TO REMOVE
			else if (zij_cij == 0.0f && Max == 0.0f && AlreadyCalced == MCF_FALSE)
			{
				Max = zij_cij;
				M = l;
				AlreadyCalced = MCF_TRUE;
			}
			else if (zij_cij == Max && AlreadyCalced == MCF_TRUE)
			{
				//choose the one with the smallest unit cost
				if (Cij_l[l] < Cij_l[M])
				{
					Max = zij_cij;
					M = l;
					AlreadyCalced = MCF_TRUE;
				}
			}
			//MIGHT NEED TO REMOVE THE ABOVE
			if (zij_cij < 0)
			{
				NegCount++;
			}
		}
	}

	printf("\nChoosing Link %i, with Value %f", M, Max);

	if (NegCount == BasisCount)
	{
		ConvLinkCandidate CLC =
		{
			.Converged = 1,
			.LinkCandidate = M,
			.CandidateValue = Max
		};
		return &CLC;
	}
	else
	{
		ConvLinkCandidate CLC =
		{
			.Converged = 0,
			.LinkCandidate = M,
			.CandidateValue = Max
		};
		return &CLC;
	}
}

//flows also need to be positive
MCFBool balanced_Network(const unsigned int N, const unsigned int L, const float* B_n, const MCFBool* Basis_l, const unsigned int* Ni_l, const unsigned int* Nj_l, const float* Xij_l)
{
	//check to make sure we have all positive values for the flows or else they will be reversing traffic!!
	MCFBool AllPos = MCF_TRUE;
	for (unsigned int l = 0; l < L; l++)
	{
		//must check IF BASIS
		if (Basis_l[l] == MCF_TRUE)
		{
			//soft constraint again
			float l_val = Xij_l[1];
			if (l_val < -0.1f)
				AllPos = MCF_FALSE;
		}
	}

	if (AllPos == MCF_FALSE)
		return MCF_FALSE;
	else
	{
		float* current_B_n = malloc(N * sizeof(float));
		for (unsigned int n = 0; n < N; n++)
		{
			current_B_n[n] = 0.0f;
		}

		for (unsigned int l = 0; l < L; l++)
		{
			if (Basis_l[l] == MCF_TRUE)
			{
				float xij_l = Xij_l[l];

				unsigned int ni_l = Ni_l[l];
				unsigned int nj_l = Nj_l[l];
				//add flow to the node we are leaving
				//subtract flow from the node we are entering
				if (nj_l > ni_l)
				{
					current_B_n[ni_l] += xij_l;
					current_B_n[nj_l] -= xij_l;
				}
				else
				{
					current_B_n[ni_l] -= xij_l;
					current_B_n[nj_l] += xij_l;
				}
			}
		}

		MCFBool Satisfied = MCF_TRUE;
		for (unsigned int n = 0; n < N; n++)
		{
			//put in a soft constrait
			if (current_B_n[n] < (B_n[n] - 0.1f) || current_B_n[n] > (B_n[n] + 0.1f))
				Satisfied = MCF_FALSE;
		}

		free(current_B_n);
		return Satisfied;
	}
}

void addBasis(const unsigned int L, const unsigned int LinkCandidate, MCFBool* Basis_1)
{
	printf("\n\nTrying to add a basis:");
	for (unsigned int l = 0; l < L; l++)
	{
		if (LinkCandidate == l)
		{
			MCFBool Was = Basis_1[l];
			Basis_1[l] = MCF_TRUE;
			printf("\nBasis %i bool was %i, and now is %i", l, Was, Basis_1[l]);
		}

	}
}

MCFBool* alloc_DeadEndNodes(const unsigned int N, const unsigned int L, const unsigned int* Ni_l, const unsigned int* Nj_l, const MCFBool* Basis_l)
{
	unsigned int* LinkCounts_n = malloc(N * sizeof(unsigned int));
	for (unsigned int n = 0; n < N; n++)
	{
		LinkCounts_n[n] = 0;
	}

	for (unsigned int l = 0; l < L; l++)
	{
		//gotta make sure we are looking at the basis variables
		if (Basis_l[l] == MCF_TRUE)
		{
			unsigned int ni_l = Ni_l[l];
			unsigned int nj_l = Nj_l[l];
			LinkCounts_n[ni_l] += 1;
			LinkCounts_n[nj_l] += 1;
		}
	}

	MCFBool* DeadEnd_n = malloc(N * sizeof(MCFBool));
	//printf("\n\nPrinting information about dead end nodes");
	for (unsigned int n = 0; n < N; n++)
	{
		if (LinkCounts_n[n] >= 2)
			DeadEnd_n[n] = MCF_FALSE;
		else
			DeadEnd_n[n] = MCF_TRUE;
		//printf("\nNode %i, Dead End Node bool: %i with count: %i", n, DeadEnd_n[n], LinkCounts_n[n]);
	}

	free(LinkCounts_n);
	return DeadEnd_n;
}

unsigned int* alloc_ActiveLinkCounts(const unsigned int N, const unsigned int L, const unsigned int* Ni_l, const unsigned int* Nj_l, const MCFBool* Basis_l)
{
	unsigned int* Counts = malloc(N * sizeof(unsigned));
	for (unsigned int n = 0; n < N; n++)
		Counts[n] = 0;

	for (unsigned int l = 0; l < L; l++)
	{
		if (Basis_l[l] == MCF_TRUE)
		{
			unsigned int ni_l = Ni_l[l];
			unsigned int nj_l = Nj_l[l];
			Counts[ni_l] += 1;
			Counts[nj_l] += 1;
		}
	}
	return Counts;
}

//static meaning they won't update during the solving phase
MCFBool* alloc_StaticLinks(const unsigned int N, const unsigned int L, const unsigned int* Ni_l, const unsigned int* Nj_l, const MCFBool* DeadEnd_n, const MCFBool* Basis_l)
{
	//printf("\n\nPrinting information about static (with a dead end) links");
	MCFBool* StaticLinks_l = malloc(L * sizeof(MCFBool));
	for (unsigned int l = 0; l < L; l++)
	{
		unsigned int ni_l = Ni_l[l];
		unsigned int nj_l = Nj_l[l];
		MCFBool Dead_i = DeadEnd_n[ni_l];
		MCFBool Dead_j = DeadEnd_n[nj_l];

		if (Dead_i == MCF_TRUE || Dead_j == MCF_TRUE)
			StaticLinks_l[l] = MCF_TRUE;
		else
			StaticLinks_l[l] = MCF_FALSE;

		//printf("\nLink %i, with ni_l %i, and nj_l %i, and static link bool: %i", l, ni_l, nj_l, StaticLinks_l[l]);
	}

	unsigned int* LinkCounts = alloc_ActiveLinkCounts(N, L, Ni_l, Nj_l, Basis_l);
	//WARNING I ADDED THIS TO CODE THAT WORKED IN SOME CASES ALREADY HOPEFULLY IT SURVIVES THIS
	//REMOVE IF NEEDED
	//need to find a not yet static link that shares a node with a static link AND THAT NODE ONLY HAS 2 LINKS
	for (unsigned int l = 0; l < L; l++)
	{
		if (StaticLinks_l[l] == MCF_FALSE && Basis_l[l] == MCF_TRUE)
		{
			unsigned int ni_l = Ni_l[l];
			unsigned int nj_l = Nj_l[l];
			for (unsigned int m = 0; m < L; m++)
			{
				if (StaticLinks_l[m] == MCF_TRUE && Basis_l[m] == MCF_TRUE && m != l)
				{
					unsigned int ni_m = Ni_l[m];
					unsigned int nj_m = Nj_l[m];
					if (ni_l == ni_m || ni_l == nj_m)
					{
						unsigned int LCount = LinkCounts[ni_l];
						if (LCount < 3)
							StaticLinks_l[l] = MCF_TRUE;
					}

					if (nj_l == ni_m || nj_l == nj_m)
					{
						unsigned int LCount = LinkCounts[nj_l];
						if (LCount < 3)
							StaticLinks_l[l] = MCF_TRUE;
					}

				}
			}
		}
	}

	free(LinkCounts);

	return StaticLinks_l;
}

void reset_LinkCalc(const unsigned int L, const MCFBool* StaticLinks_l, const float* Xij_l, float* LinkFlow_l, MCFBool* CalcedL_l, const unsigned int Theta_l, const float Theta)
{
	//printf("\n\nPrinting information on whether or not Link is precalculated because of static links and theta");
	for (unsigned int l = 0; l < L; l++)
	{
		if (StaticLinks_l[l] == MCF_TRUE)
		{
			LinkFlow_l[l] = Xij_l[l];
			CalcedL_l[l] = MCF_TRUE;
		}
		else if (l == Theta_l)
		{
			LinkFlow_l[l] = Theta;
			CalcedL_l[l] = MCF_TRUE;
		}
		else
		{
			LinkFlow_l[l] = 0.0f;
			CalcedL_l[l] = MCF_FALSE;
		}
		//printf("\nLink %i, from %i to %i has calced static bool %i with value %f", l CalcedL_l[l], LinkFlow_l[l]);
	}
}

MCFBool doneCalculating_l(const unsigned int L, const MCFBool* CalcedL_l, const MCFBool* Basis_l)
{
	MCFBool Done_l = MCF_TRUE;
	for (unsigned int l = 0; l < L; l++)
	{
		if (Basis_l[l] == MCF_TRUE)
		{
			if (CalcedL_l[l] == MCF_FALSE)
			{
				Done_l = MCF_FALSE;
			}

		}
	}
	//if (Done_l == MCF_FALSE)
	//	printf("\nNot done calculating links\n");
	//else
	//	printf("\nFinished calculating links\n");
	return Done_l;
}


float* try_Theta(const unsigned int N, const unsigned int L, const unsigned int* Ni_l, const unsigned int* Nj_l, const MCFBool* StaticLinks_l,
	const float* B_n, const float* Xij_l, const float Theta, const int Theta_l, const MCFBool* Basis_l)
{
	float* LinkFlow_l = malloc(L * sizeof(float));
	MCFBool* CalcedL_l = malloc(L * sizeof(MCFBool));
	reset_LinkCalc(L, StaticLinks_l, Xij_l, LinkFlow_l, CalcedL_l, Theta_l, Theta);

	//While we are not done.  DONE = True if(ALL BASIS LINKS CALCULATED && all STATIC BASIS should already be considered CALCULATED)
	while (doneCalculating_l(L, CalcedL_l, Basis_l) == MCF_FALSE)
	{
		//LOOP for each link that IS a BASIS LINK && NOT CALCULATED && NOT STATIC
		for (unsigned int l = 0; l < L; l++)
		{
			if (Basis_l[l] == MCF_TRUE && CalcedL_l[l] == MCF_FALSE && StaticLinks_l[l] == MCF_FALSE)
			{
				//for each of these find if we touch a: CALCULATED, NON-STATIC, BASIS
				//if we find multiple we will just overwrite with las
				//have to loop again because I don't have a map from n to l (neighbor l assigned to m.  make sure it's a good neighbor!!!
				unsigned int ni_l = Ni_l[l];
				unsigned int nj_l = Nj_l[l];

				//neighbor information
				MCFBool GoodNeighbor = MCF_FALSE;
				unsigned int n_m = N; //shared node
				unsigned int l_m = N; //link index of match
				float xn_m = 0; //value relative to the node that we are considering

				for (unsigned int m = 0; m < L; m++)
				{
					if (GoodNeighbor == MCF_FALSE) //repetition should work but it is causing some nasty bugs to appear in my program
					{
						if (CalcedL_l[m] == MCF_TRUE && StaticLinks_l[m] == MCF_FALSE && Basis_l[m] == MCF_TRUE)
						{
							unsigned int ni_m = Ni_l[m];
							unsigned int nj_m = Nj_l[m];
							float xi_m = 0;
							float xj_m = 0;
							//if leaving i set positive for i and negative for j. else negative for i and positive for j
							if (ni_m < nj_m)
							{
								xi_m = LinkFlow_l[m];
								xj_m = -LinkFlow_l[m];
							}
							else
							{
								xi_m = -LinkFlow_l[m];
								xj_m = LinkFlow_l[m];
							}

							//more if
							if (ni_l == ni_m)
							{
								GoodNeighbor = MCF_TRUE;
								l_m = m;
								n_m = ni_l;  //==ni_Candidate
								xn_m = xi_m;
							}
							else if (ni_l == nj_m)
							{
								GoodNeighbor = MCF_TRUE;
								l_m = m;
								n_m = ni_l;  //==nj_Candidate
								xn_m = xj_m;
							}
							if (nj_l == ni_m)
							{
								GoodNeighbor = MCF_TRUE;
								l_m = m;
								n_m = nj_l; //==ni_Candidate
								xn_m = xi_m;
							}
							else if (nj_l == nj_m)
							{
								GoodNeighbor = MCF_TRUE;
								l_m = m;
								n_m = nj_l;  //==nj_Candidate
								xn_m = xj_m;
							}
						}
					}
				}
				//alright that's done...
				if (GoodNeighbor == MCF_TRUE)
				{
					float n_Val = B_n[n_m];
					float left_Over = n_Val - xn_m;

					MCFBool Negative = MCF_FALSE;
					if (n_m == ni_l && ni_l > nj_l)
						Negative = MCF_TRUE;
					if (n_m == nj_l && nj_l > ni_l)
						Negative = MCF_TRUE;

					LinkFlow_l[l] = abs(left_Over);
					CalcedL_l[l] = MCF_TRUE;

					////usiing soft constraints
					//if (left_Over < -0.1f && Negative == MCF_TRUE)
					//{
					//	LinkFlow_l[l] = abs(left_Over);
					//	CalcedL_l[l] = MCF_TRUE;
					//	//printf("\nValue of link %i is set to %f\n", l, LinkFlow_l[l]);
					//}
					//else if (left_Over > 0.1f && Negative == MCF_FALSE)
					//{
					//	LinkFlow_l[l] = left_Over;
					//	CalcedL_l[l] = MCF_TRUE;
					//	//printf("\nValue of link %i is set to %f\n", l, LinkFlow_l[l]);
					//}
					//else if (left_Over <= 0.1f && left_Over >= -0.1f)
					//{
					//	LinkFlow_l[l] = 0.0f;
					//	CalcedL_l[l] = MCF_TRUE;
					//	//printf("\nValue of link %i is set to %f \n", l, LinkFlow_l[l]);
					//}

					////added this
					//else if (left_Over > 0.1f && Negative == MCF_FALSE)
					//{
					//	LinkFlow_l[l] = left_Over;
					//	CalcedL_l[l] = MCF_TRUE;
					//	//printf("\nValue of link %i is set to %f\n", l, LinkFlow_l[l]);
					//}

					//else
					//{
					//	printf("\n\nSomething went wrong with the value of the newly calculated link AND PRINTING DEBUG CALCULATION DATA\n");
					//	for (unsigned int l = 0; l < L; l++)
					//	{
					//		printf("\nLink %i, ni_l %i, nj_l %i, static bool: %i", l, Ni_l[l], Nj_l[l], StaticLinks_l[l]);
					//	}
					//	getchar();
					//}
				}
			}
		}
	}
	free(CalcedL_l);
	return(LinkFlow_l);
	//end of loop
}

unsigned int count_True(unsigned int L, MCFBool* True_l)
{
	unsigned int count = 0;
	for (unsigned int l = 0; l < L; l++)
	{
		if (True_l[l] == MCF_TRUE)
			count++;
	}
	return count;
}

unsigned int count_NewBasis(const unsigned int L, const MCFBool* basis_l, const float* LinkFlow_l)
{
	unsigned int count = 0;
	for (unsigned int l = 0; l < L; l++)
	{
		if (basis_l[l] == MCF_TRUE)
		{
			if (LinkFlow_l[l] > 0.1f || LinkFlow_l[l] < -0.1f)
				count++;
		}
	}
	return count;
}

float* alloc_B_n_OnStatic(unsigned int N, unsigned int L, unsigned int* Ni_l, unsigned int* Nj_l, MCFBool* StaticLinks_l, const float* B_n, const float* Xij_l)
{
	float* dist_B_n = malloc(N * sizeof(float));
	for (unsigned int n = 0; n < N; n++)
	{
		dist_B_n[n] = B_n[n];
	}

	//the new node values should be what they now desire
	//calculate this as B_n - known
	for (unsigned int l = 0; l < L; l++)
	{
		//only for knows (i.e. static links)
		if (StaticLinks_l[l] == MCF_TRUE)
		{
			unsigned int ni_l = Ni_l[l];
			unsigned int nj_l = Nj_l[l];
			float vali_l = 0;
			float valj_l = 0;
			if (ni_l < nj_l)
			{
				vali_l = Xij_l[l];
				valj_l = -Xij_l[l];
			}

			else
			{
				vali_l = -Xij_l[l];
				valj_l = Xij_l[l];
			}
			dist_B_n[ni_l] = B_n[ni_l] - vali_l;
			dist_B_n[nj_l] = B_n[nj_l] - valj_l;
		}
	}

	//printf("\nPrinting redistributed B_n values based on static links:");
	//for (unsigned int n = 0; n < N; n++)
	//{
	//	printf("\nNodes %i value is now: %f", n, dist_B_n[n]);
	//}

	return dist_B_n;
}

void calcNB(const unsigned int N, const unsigned int L, const unsigned int* Ni_l, const unsigned int* Nj_l, const float* B_n, MCFBool* Basis_l, float* Xij_l, const unsigned int LinkCandidate)
{
	addBasis(L, LinkCandidate, Basis_l);
	//find the dead end nodes
	MCFBool* DeadEnd_n = alloc_DeadEndNodes(N, L, Ni_l, Nj_l, Basis_l); //also 
	MCFBool* StaticLinks_l = alloc_StaticLinks(N, L, Ni_l, Nj_l, DeadEnd_n, Basis_l);
	MCFBool* dist_B_n = alloc_B_n_OnStatic(N, L, Ni_l, Nj_l, StaticLinks_l, B_n, Xij_l);

	unsigned int CountBasis = count_True(L, Basis_l);
	unsigned int CountNewBasis = CountBasis;
	unsigned int ThetaInt = 0;
	while (CountNewBasis >= CountBasis)
	{
		ThetaInt++;
		float Theta = 1.0f * ThetaInt;
		unsigned int Theta_l = LinkCandidate;

		float* LinkFlow_l = try_Theta(N, L, Ni_l, Nj_l, StaticLinks_l, dist_B_n, Xij_l, Theta, Theta_l, Basis_l);
		CountNewBasis = count_NewBasis(L, Basis_l, LinkFlow_l);

		//printf("\nNumber of new Bases: %i\n", CountNewBasis);
		free(LinkFlow_l);
	}

	float Theta = 1.0f * ThetaInt;
	unsigned int Theta_l = LinkCandidate;
	float* LinkFlow_l = try_Theta(N, L, Ni_l, Nj_l, StaticLinks_l, dist_B_n, Xij_l, Theta, Theta_l, Basis_l);

	printf("\n\nThe found value of Theta is %f\n", Theta);

	//assign values from LinkFlow to Xij
	for (unsigned int l = 0; l < L; l++)
	{
		Xij_l[l] = LinkFlow_l[l];
		if (LinkFlow_l[l] > 0.1f || LinkFlow_l[l] < -0.1f)
			Basis_l[l] = MCF_TRUE;
		else
			Basis_l[l] = MCF_FALSE;
	}

	unsigned int NumNewBasis = count_True(L, Basis_l);
	unsigned int NumActiveLinks = count_NewBasis(L, Basis_l, LinkFlow_l);
	if (NumNewBasis == NumActiveLinks)
		printf("\nThe number of non-zero links and the number of new bases match after calcNB with count: %i\n", NumNewBasis);
	else
		printf("\nThe number of non-zero links is %i and the number of new basis is %i.  THESE DON'T MATCH.\n", NumActiveLinks, NumNewBasis);

	printf("\nThe current values for the links:");
	for (unsigned int l = 0; l < L; l++)
	{
		printf("\nLink %i, from node %i to node %i with value: %f and basis bool: %i", l, Ni_l[l], Nj_l[l], Xij_l[l], Basis_l[l]);
	}

	free(LinkFlow_l);
	free(DeadEnd_n);
	free(StaticLinks_l);
	free(dist_B_n);
}

void phaseII_Iterate2(const unsigned int N, const unsigned int L, const unsigned int* Ni_l, const unsigned int* Nj_l, const float* Cij_l, const float* B_n, MCFBool* Basis_l, float* Xij_l)
{

	float currentCost = cost(L, Cij_l, Xij_l);
	unsigned int Iteration = 0;
	printf("\n\nThe current cost at iteration 0 is: %f", currentCost);
	for (unsigned iter = 0; iter < 10; iter++)
	{
		float* W_n = compW(N, L, Basis_l, Ni_l, Nj_l, Cij_l);

		ConvLinkCandidate* CLC = create_ConvLinkCandidate(N, L, Basis_l, Ni_l, Nj_l, W_n, Cij_l);
		unsigned int LinkCandidate = CLC->LinkCandidate;
		calcNB(N, L, Ni_l, Nj_l, B_n, Basis_l, Xij_l, LinkCandidate);

		currentCost = cost(L, Cij_l, Xij_l);
		printf("\n\nThe current cost at iteration %i is: %f", Iteration, currentCost);
		free(W_n);
		getchar();
	}
}

//need starting flows and starting basis variables (same links)
void phaseII_Iterate(const unsigned int N, const unsigned int L, const unsigned int* Ni_l, const unsigned int* Nj_l, const float* Cij_l, const float* B_n, MCFBool* Basis_l, float* Xij_l)
{
	for (unsigned int i = 0; i < 100; i++)
	{
		printf("\n\nIteration %i", i);
		float currentCost = cost(L, Cij_l, Xij_l);
		printf("\nThe current cost is: %f", currentCost);
		float* W_n = compW(N, L, Basis_l, Ni_l, Nj_l, Cij_l);

		ConvLinkCandidate* CLC = create_ConvLinkCandidate(N, L, Basis_l, Ni_l, Nj_l, W_n, Cij_l);  //don't attempt to free a stack pointer

		free(W_n);

		MCFBool Converged = CLC->Converged;

		if (Converged != MCF_FALSE)
		{
			float finalCost = cost(L, Cij_l, Xij_l);
			printf("\n\nWe have converged after %i steps with final cost: %f\n", i, finalCost);
			printf("\nThe link values are:");
			for (unsigned int l = 0; l < L; l++)
			{
				printf("\nLink %i has value %f", l, Xij_l[l]);
			}
			getchar();
		}
		else
		{
			unsigned int LinkCandidate = CLC->LinkCandidate;
			//update_Basis(N, L, LinkCandidate, B_n, Ni_l, Nj_l, Basis_l, Xij_l);
			//
			calcNB(N, L, Ni_l, Nj_l, B_n, Basis_l, Xij_l, LinkCandidate);
		}
	}
}
//
void phaseII_problem_example()
{
	printf("\n\nPhase II problem\n");
	unsigned int N = 4;
	float* B_n = malloc(N * sizeof(float));

	B_n[0] = 4.0f;
	B_n[1] = 2.0f;
	B_n[2] = -1.0f;
	B_n[3] = -5.0f;

	//put the check for B_n sums here

	unsigned L = 7;
	MCFBool* Basis_l = malloc(L * sizeof(MCFBool));
	unsigned int* Ni_l = malloc(L * sizeof(unsigned int));
	unsigned int* Nj_l = malloc(L * sizeof(unsigned int));
	float* Cij_l = malloc(L * sizeof(float));
	float* Xij_l = malloc(L * sizeof(float));

	//don't forget we are counting in memory OFFSET!!!
	Basis_l[0] = MCF_TRUE;
	Basis_l[1] = MCF_TRUE;
	Basis_l[2] = MCF_FALSE;
	Basis_l[3] = MCF_TRUE;
	Basis_l[4] = MCF_FALSE;
	Basis_l[5] = MCF_FALSE;
	Basis_l[6] = MCF_FALSE;

	Ni_l[0] = 0;
	Ni_l[1] = 0;
	Ni_l[2] = 3;
	Ni_l[3] = 1;
	Ni_l[4] = 1;
	Ni_l[5] = 2;
	Ni_l[6] = 2;

	Nj_l[0] = 1;
	Nj_l[1] = 2;
	Nj_l[2] = 0;
	Nj_l[3] = 3;
	Nj_l[4] = 2;
	Nj_l[5] = 1;
	Nj_l[6] = 3;

	Cij_l[0] = 2.0f;
	Cij_l[1] = -5.0f;
	Cij_l[2] = 7.0f;
	Cij_l[3] = 4.0f;
	Cij_l[4] = -1.0f;
	Cij_l[5] = 6.0f;
	Cij_l[6] = 3.0f;

	Xij_l[0] = 3.0f;
	Xij_l[1] = 1.0f;
	Xij_l[2] = 0.0f;
	Xij_l[3] = 5.0f;
	Xij_l[4] = 0.0f;
	Xij_l[5] = 0.0f;
	Xij_l[6] = 0.0f;

	
	phaseII_Iterate2(N, L, Ni_l, Nj_l, Cij_l, B_n, Basis_l, Xij_l);

	free(B_n);

	free(Basis_l);
	free(Ni_l);
	free(Nj_l);
	free(Cij_l);
	free(Xij_l);
}

void phaseII_problem_3()
{
	printf("\n\nPhase II problem\n");
	unsigned int N = 8;
	float* B_n = malloc(N * sizeof(float));

	B_n[0] = 10.0f;
	B_n[1] = 20.0f;
	B_n[2] = 0.0f;
	B_n[3] = -5.0f;
	B_n[4] = 0.0f;
	B_n[5] = 0.0f;
	B_n[6] = -15.0f;
	B_n[7] = -10.0f;

	//put the check for B_n sums here

	unsigned L = 11;
	MCFBool* Basis_l = malloc(L * sizeof(MCFBool));
	unsigned int* Ni_l = malloc(L * sizeof(unsigned int));
	unsigned int* Nj_l = malloc(L * sizeof(unsigned int));
	float* Cij_l = malloc(L * sizeof(float));
	float* Xij_l = malloc(L * sizeof(float));

	//don't forget we are counting in memory OFFSET!!!
	Basis_l[0] = MCF_TRUE;
	Basis_l[1] = MCF_FALSE;
	Basis_l[2] = MCF_TRUE;
	Basis_l[3] = MCF_TRUE;
	Basis_l[4] = MCF_TRUE;
	Basis_l[5] = MCF_TRUE;
	Basis_l[6] = MCF_TRUE;
	Basis_l[7] = MCF_FALSE;
	Basis_l[8] = MCF_TRUE;
	Basis_l[9] = MCF_TRUE;
	Basis_l[10] = MCF_FALSE;

	Ni_l[0] = 0;
	Ni_l[1] = 1;
	Ni_l[2] = 1;
	Ni_l[3] = 1;
	Ni_l[4] = 2;
	Ni_l[5] = 2;
	Ni_l[6] = 3;
	Ni_l[7] = 4;
	Ni_l[8] = 4;
	Ni_l[9] = 5;
	Ni_l[10] = 6;

	Nj_l[0] = 3;
	Nj_l[1] = 0;
	Nj_l[2] = 2;
	Nj_l[3] = 5;
	Nj_l[4] = 3;
	Nj_l[5] = 4;
	Nj_l[6] = 6;
	Nj_l[7] = 5;
	Nj_l[8] = 6;
	Nj_l[9] = 7;
	Nj_l[10] = 7;

	Cij_l[0] = 2.0f;
	Cij_l[1] = 1.0f;
	Cij_l[2] = 0.0f;
	Cij_l[3] = 6.0f;
	Cij_l[4] = 1.0f;
	Cij_l[5] = 4.0f;
	Cij_l[6] = 5.0f;
	Cij_l[7] = 2.0f;
	Cij_l[8] = 7.0f;
	Cij_l[9] = 8.0f;
	Cij_l[10] = 9.0f;

	Xij_l[0] = 10.0f;
	Xij_l[1] = 0.0f;
	Xij_l[2] = 10.0f;
	Xij_l[3] = 10.0f;
	Xij_l[4] = 5.0f;
	Xij_l[5] = 5.0f;
	Xij_l[6] = 10.0f;
	Xij_l[7] = 0.0f;
	Xij_l[8] = 5.0f;
	Xij_l[9] = 10.0f;
	Xij_l[10] = 0.0f;

	phaseII_Iterate2(N, L, Ni_l, Nj_l, Cij_l, B_n, Basis_l, Xij_l);

	free(B_n);

	free(Basis_l);
	free(Ni_l);
	free(Nj_l);
	free(Cij_l);
	free(Xij_l);
}

void phaseII_problem()
{
	printf("\n\nPhase II problem\n");
	unsigned int N = 4;
	float* B_n = malloc(N * sizeof(float));

	B_n[0] = 3.0f;
	B_n[1] = 3.0f;
	B_n[2] = 1.0f;
	B_n[3] = 7.0f;

	//put the check for B_n sums here

	unsigned L = 5;
	MCFBool* Basis_l = malloc(L * sizeof(MCFBool));
	unsigned int* Ni_l = malloc(L * sizeof(unsigned int));
	unsigned int* Nj_l = malloc(L * sizeof(unsigned int));
	float* Cij_l = malloc(L * sizeof(float));
	float* Xij_l = malloc(L * sizeof(float));

	//don't forget we are counting in memory OFFSET!!!
	Basis_l[0] = MCF_TRUE;  //3
	Basis_l[1] = MCF_FALSE;  //0
	Basis_l[2] = MCF_TRUE;  //7
	Basis_l[3] = MCF_TRUE;  //1
	Basis_l[4] = MCF_FALSE;  //0

	Ni_l[0] = 0;
	Ni_l[1] = 0;
	Ni_l[2] = 1;
	Ni_l[3] = 2;
	Ni_l[4] = 2;

	Nj_l[0] = 1;
	Nj_l[1] = 2;
	Nj_l[2] = 3;
	Nj_l[3] = 1;
	Nj_l[4] = 3;

	Cij_l[0] = 3.0f;
	Cij_l[1] = 4.0f;
	Cij_l[2] = 2.0f;
	Cij_l[3] = -1.0f;
	Cij_l[4] = 6.0f;

	Xij_l[0] = 3.0f;
	Xij_l[1] = 0.0f;
	Xij_l[2] = 7.0f;
	Xij_l[3] = 1.0f;
	Xij_l[4] = 0.0f;

	phaseII_Iterate(N, L, Ni_l, Nj_l, Cij_l, B_n, Basis_l, Xij_l);

	free(B_n);

	free(Basis_l);
	free(Ni_l);
	free(Nj_l);
	free(Cij_l);
	free(Xij_l);
}

void phaseII_problem_4()
{
	printf("\n\nPhase II problem\n");
	unsigned int N = 33;
	float* B_n = malloc(N * sizeof(float));

	B_n[0] = 20.0f;
	B_n[1] = 0.0f;
	B_n[2] = 0.0f;
	B_n[3] = 0.0f;
	B_n[4] = 0.0f;
	B_n[5] = 0.0f;
	B_n[6] = 0.0f;
	B_n[7] = 0.0f;
	B_n[8] = 0.0f;
	B_n[9] = 0.0f;
	B_n[10] = 0.0f;
	B_n[11] = 0.0f;
	B_n[12] = 0.0f;
	B_n[13] = 0.0f;
	B_n[14] = 0.0f;
	B_n[15] = 0.0f;
	B_n[16] = -5.0f;
	B_n[17] = 0.0f;
	B_n[18] = 0.0f;
	B_n[19] = 0.0f;
	B_n[20] = -5.0f;
	B_n[21] = -5.0f;
	B_n[22] = 0.0f;
	B_n[23] = 0.0f;
	B_n[24] = 0.0f;
	B_n[25] = 0.0f;
	B_n[26] = 0.0f;
	B_n[27] = 0.0f;
	B_n[28] = 0.0f;
	B_n[29] = 0.0f;
	B_n[30] = 0.0f;
	B_n[31] = 0.0f;
	B_n[32] = -5.0f;

	//put the check for B_n sums here

	unsigned L = 51;
	MCFBool* Basis_l = malloc(L * sizeof(MCFBool));
	unsigned int* Ni_l = malloc(L * sizeof(unsigned int));
	unsigned int* Nj_l = malloc(L * sizeof(unsigned int));
	float* Cij_l = malloc(L * sizeof(float));
	float* Xij_l = malloc(L * sizeof(float));

	//don't forget we are counting in memory OFFSET!!!
	Basis_l[0] = MCF_FALSE;
	Basis_l[1] = MCF_TRUE;
	Basis_l[2] = MCF_FALSE;
	Basis_l[3] = MCF_FALSE;
	Basis_l[4] = MCF_TRUE;
	Basis_l[5] = MCF_FALSE;
	Basis_l[6] = MCF_FALSE;
	Basis_l[7] = MCF_TRUE;
	Basis_l[8] = MCF_FALSE;
	Basis_l[9] = MCF_FALSE;

	Basis_l[10] = MCF_TRUE;
	Basis_l[11] = MCF_FALSE;
	Basis_l[12] = MCF_FALSE;
	Basis_l[13] = MCF_FALSE;
	Basis_l[14] = MCF_FALSE;
	Basis_l[15] = MCF_FALSE;
	Basis_l[16] = MCF_FALSE;
	Basis_l[17] = MCF_FALSE;
	Basis_l[18] = MCF_TRUE;
	Basis_l[19] = MCF_TRUE;

	Basis_l[20] = MCF_FALSE;
	Basis_l[21] = MCF_TRUE;
	Basis_l[22] = MCF_FALSE;
	Basis_l[23] = MCF_FALSE;
	Basis_l[24] = MCF_TRUE;
	Basis_l[25] = MCF_FALSE;
	Basis_l[26] = MCF_FALSE;
	Basis_l[27] = MCF_FALSE;
	Basis_l[28] = MCF_FALSE;
	Basis_l[29] = MCF_FALSE;

	Basis_l[30] = MCF_FALSE;
	Basis_l[31] = MCF_FALSE;
	Basis_l[32] = MCF_FALSE;
	Basis_l[33] = MCF_FALSE;
	Basis_l[34] = MCF_FALSE;
	Basis_l[35] = MCF_FALSE;
	Basis_l[36] = MCF_FALSE;
	Basis_l[37] = MCF_FALSE;
	Basis_l[38] = MCF_FALSE;
	Basis_l[39] = MCF_FALSE;

	Basis_l[40] = MCF_FALSE;
	Basis_l[41] = MCF_FALSE;
	Basis_l[42] = MCF_FALSE;
	Basis_l[43] = MCF_FALSE;
	Basis_l[44] = MCF_FALSE;
	Basis_l[45] = MCF_FALSE;
	Basis_l[46] = MCF_FALSE;
	Basis_l[47] = MCF_FALSE;
	Basis_l[48] = MCF_FALSE;
	Basis_l[49] = MCF_FALSE;

	Basis_l[50] = MCF_FALSE;

	Ni_l[0] = 0;
	Ni_l[1] = 0;
	Ni_l[2] = 2;
	Ni_l[3] = 2;
	Ni_l[4] = 4;
	Ni_l[5] = 4;
	Ni_l[6] = 4;
	Ni_l[7] = 6;
	Ni_l[8] = 5;
	Ni_l[9] = 18;

	Nj_l[0] = 1;
	Nj_l[1] = 2;
	Nj_l[2] = 1;
	Nj_l[3] = 3;
	Nj_l[4] = 6;
	Nj_l[5] = 3;
	Nj_l[6] = 5;
	Nj_l[7] = 7;
	Nj_l[8] = 10;
	Nj_l[9] = 17;

	Cij_l[0] = 4;
	Cij_l[1] = 5;
	Cij_l[2] = 2;
	Cij_l[3] = 7;
	Cij_l[4] = 8;
	Cij_l[5] = 9;
	Cij_l[6] = 5;
	Cij_l[7] = 2;
	Cij_l[8] = 7;
	Cij_l[9] = 8;

	Xij_l[0] = 0;
	Xij_l[1] = 20;
	Xij_l[2] = 0;
	Xij_l[3] = 0;
	Xij_l[4] = 20;
	Xij_l[5] = 0;
	Xij_l[6] = 0;
	Xij_l[7] = 20;
	Xij_l[8] = 0;
	Xij_l[9] = 0;

	Ni_l[10] = 7;
	Ni_l[11] = 6;
	Ni_l[12] = 5;
	Ni_l[13] = 9;
	Ni_l[14] = 6;
	Ni_l[15] = 8;
	Ni_l[16] = 9;
	Ni_l[17] = 11;
	Ni_l[18] = 18;
	Ni_l[19] = 20;

	Nj_l[10] = 20;
	Nj_l[11] = 19;
	Nj_l[12] = 8;
	Nj_l[13] = 10;
	Nj_l[14] = 7;
	Nj_l[15] = 11;
	Nj_l[16] = 10;
	Nj_l[17] = 12;
	Nj_l[18] = 16;
	Nj_l[19] = 18;

	Cij_l[10] = 9;
	Cij_l[11] = 2;
	Cij_l[12] = 1;
	Cij_l[13] = 0;
	Cij_l[14] = 6;
	Cij_l[15] = 1;
	Cij_l[16] = 7;
	Cij_l[17] = 4;
	Cij_l[18] = 4;
	Cij_l[19] = 7;

	Xij_l[10] = 20;
	Xij_l[11] = 0;
	Xij_l[12] = 0;
	Xij_l[13] = 0;
	Xij_l[14] = 0;
	Xij_l[15] = 0;
	Xij_l[16] = 0;
	Xij_l[17] = 0;
	Xij_l[18] = 5;
	Xij_l[19] = 5;

	Ni_l[20] = 13;
	Ni_l[21] = 20;
	Ni_l[22] = 12;
	Ni_l[23] = 13;
	Ni_l[24] = 21;
	Ni_l[25] = 17;
	Ni_l[26] = 14;
	Ni_l[27] = 15;
	Ni_l[28] = 21;
	Ni_l[29] = 17;

	Nj_l[20] = 14;
	Nj_l[21] = 21;
	Nj_l[22] = 15;
	Nj_l[23] = 16;
	Nj_l[24] = 32;
	Nj_l[25] = 18;
	Nj_l[26] = 19;
	Nj_l[27] = 20;
	Nj_l[28] = 16;
	Nj_l[29] = 18;

	Cij_l[20] = 14;
	Cij_l[21] = 15;
	Cij_l[22] = 10;
	Cij_l[23] = 7;
	Cij_l[24] = 5;
	Cij_l[25] = 6;
	Cij_l[26] = 7;
	Cij_l[27] = 8;
	Cij_l[28] = 10;
	Cij_l[29] = 9;

	Xij_l[20] = 0;
	Xij_l[21] = 10;
	Xij_l[22] = 0;
	Xij_l[23] = 0;
	Xij_l[24] = 5;
	Xij_l[25] = 0;
	Xij_l[26] = 0;
	Xij_l[27] = 0;
	Xij_l[28] = 0;
	Xij_l[29] = 0;

	Ni_l[30] = 18;
	Ni_l[31] = 18;
	Ni_l[32] = 24;
	Ni_l[33] = 20;
	Ni_l[34] = 21;
	Ni_l[35] = 21;
	Ni_l[36] = 28;
	Ni_l[37] = 29;
	Ni_l[38] = 30;
	Ni_l[39] = 24;

	Nj_l[30] = 22;
	Nj_l[31] = 23;
	Nj_l[32] = 19;
	Nj_l[33] = 25;
	Nj_l[34] = 26;
	Nj_l[35] = 27;
	Nj_l[36] = 22;
	Nj_l[37] = 23;
	Nj_l[38] = 31;
	Nj_l[39] = 31;

	Cij_l[30] = 8;
	Cij_l[31] = 7;
	Cij_l[32] = 6;
	Cij_l[33] = 5;
	Cij_l[34] = 4;
	Cij_l[35] = 3;
	Cij_l[36] = 2;
	Cij_l[37] = 1;
	Cij_l[38] = 8;
	Cij_l[39] = 7;

	Xij_l[30] = 0;
	Xij_l[31] = 0;
	Xij_l[32] = 0;
	Xij_l[33] = 0;
	Xij_l[34] = 0;
	Xij_l[35] = 0;
	Xij_l[36] = 0;
	Xij_l[37] = 0;
	Xij_l[38] = 0;
	Xij_l[39] = 0;

	Ni_l[40] = 32;
	Ni_l[41] = 31;
	Ni_l[42] = 25;
	Ni_l[43] = 40;
	Ni_l[44] = 41;
	Ni_l[45] = 43;
	Ni_l[46] = 28;
	Ni_l[47] = 42;
	Ni_l[48] = 44;
	Ni_l[49] = 31;

	Nj_l[40] = 31;
	Nj_l[41] = 24;
	Nj_l[42] = 29;
	Nj_l[43] = 26;
	Nj_l[44] = 27;
	Nj_l[45] = 27;
	Nj_l[46] = 44;
	Nj_l[47] = 29;
	Nj_l[48] = 30;
	Nj_l[49] = 45;

	Cij_l[40] = 8;
	Cij_l[41] = 4;
	Cij_l[42] = 5;
	Cij_l[43] = 6;
	Cij_l[44] = 3;
	Cij_l[45] = 7;
	Cij_l[46] = 5;
	Cij_l[47] = 9;
	Cij_l[48] = 2;
	Cij_l[49] = 18;

	Xij_l[40] = 0;
	Xij_l[41] = 0;
	Xij_l[42] = 0;
	Xij_l[43] = 0;
	Xij_l[44] = 0;
	Xij_l[45] = 0;
	Xij_l[46] = 0;
	Xij_l[47] = 0;
	Xij_l[48] = 0;
	Xij_l[49] = 0;

	Ni_l[50] = 1;
	Nj_l[50] = 32;
	Cij_l[50] = 8;
	Xij_l[50] = 0;

	phaseII_Iterate2(N, L, Ni_l, Nj_l, Cij_l, B_n, Basis_l, Xij_l);

	free(B_n);

	free(Basis_l);
	free(Ni_l);
	free(Nj_l);
	free(Cij_l);
	free(Xij_l);
}

//stack pointer
void phaseI_Iterate(const unsigned int artN, const unsigned int artL, const unsigned int* artNi_l, const unsigned int* artNj_l, const MCFBool* artBasis_l, const float* artCij_l, const float* artXij_l, const float* artB_n)
{

	float currentCost = cost(artL, artCij_l, artXij_l);
	unsigned int Iteration = 0;
	printf("\n\nThe current cost at iteration 0 is: %f", currentCost);
	while (currentCost > 0)
	{
		Iteration++;
		float* W_n = compW(artN, artL, artBasis_l, artNi_l, artNj_l, artCij_l);

		ConvLinkCandidate* CLC = create_ConvLinkCandidate(artN, artL, artBasis_l, artNi_l, artNj_l, W_n, artCij_l);
		unsigned int LinkCandidate = CLC->LinkCandidate;
		calcNB(artN, artL, artNi_l, artNj_l, artB_n, artBasis_l, artXij_l, LinkCandidate);

		currentCost = cost(artL, artCij_l, artXij_l);
		printf("\n\nThe current cost at iteration %i is: %f", Iteration, currentCost);
		getchar();
	}
}

void phaseI_setup(unsigned int N, unsigned int L, const float* B_n, const unsigned int* Ni_l, const unsigned int* Nj_l)
{
	printf("\n\nPhase I problem\n");
	unsigned int artN = N + 1;
	float* artB_n = malloc(artN * sizeof(float));
	for (unsigned int n = 0; n < artN; n++)
	{
		if (n == N)  //the last artificial one
			artB_n[n] = 0.0f;
		else
			artB_n[n] = B_n[n];
	}

	unsigned int artL = L + N;
	MCFBool* artBasis_l = malloc(artL * sizeof(MCFBool));
	unsigned int* artNi_l = malloc(artL * sizeof(unsigned int));
	unsigned int* artNj_l = malloc(artL * sizeof(unsigned int));
	float* artCij_l = malloc(artL * sizeof(float));
	float* artXij_l = malloc(artL * sizeof(float));

	for (unsigned int l = 0; l < L; l++)
	{
		artBasis_l[l] = MCF_FALSE;
		artNi_l[l] = Ni_l[l];
		artNj_l[l] = Nj_l[l];
		artCij_l[l] = 0.0f;
		artXij_l[l] = 0.0f;
	}
	for (unsigned int n = L; n < artL; n++)
	{
		artBasis_l[n] = MCF_TRUE;
		artCij_l[n] = 1.0f;
		//assign direction based on if B_n is greater than or less than 0.0f
		if (B_n[n - L] > 0)
		{
			artNi_l[n] = n - L;
			artNj_l[n] = artN - 1;
			artXij_l[n] = B_n[n - L];
		}
		else
		{
			artNi_l[n] = artN - 1;
			artNj_l[n] = n - L;
			artXij_l[n] = -B_n[n - L];
		}
	}

	printf("\nPrinting link data for phase I setup alteration:");
	for (unsigned int l = 0; l < artL; l++)
	{
		printf("\nLink index: %i Basis Bool: %i; Node i index: %i; Node j index: %i; Link cost: %f Link flow: %f", l, artBasis_l[l], artNi_l[l], artNj_l[l], artCij_l[l], artXij_l[l]);
	}
	printf("\n\nPrinting node data of flow for phase I setup");
	for (unsigned int n = 0; n < artN; n++)
	{
		printf("\nNode %i flow value: %f", n, artB_n[n]);
	}

	phaseI_Iterate(artN, artL, artNi_l, artNj_l, artBasis_l, artCij_l, artXij_l, artB_n);

	free(artB_n);
	free(artBasis_l);
	free(artNi_l);
	free(artNj_l);
	free(artCij_l);
}

void phaseI_Problem_Example()
{
	unsigned int N = 4;
	float* B_n = malloc(N * sizeof(float));

	B_n[0] = 4.0f;
	B_n[1] = 2.0f;
	B_n[2] = -1.0f;
	B_n[3] = -5.0f;

	//put the check for B_n sums here

	unsigned L = 7;
	MCFBool* Basis_l = malloc(L * sizeof(MCFBool));
	unsigned int* Ni_l = malloc(L * sizeof(unsigned int));
	unsigned int* Nj_l = malloc(L * sizeof(unsigned int));
	float* Cij_l = malloc(L * sizeof(float));
	float* Xij_l = malloc(L * sizeof(float));

	Ni_l[0] = 0;
	Ni_l[1] = 0;
	Ni_l[2] = 3;
	Ni_l[3] = 1;
	Ni_l[4] = 1;
	Ni_l[5] = 2;
	Ni_l[6] = 2;

	Nj_l[0] = 1;
	Nj_l[1] = 2;
	Nj_l[2] = 0;
	Nj_l[3] = 3;
	Nj_l[4] = 2;
	Nj_l[5] = 1;
	Nj_l[6] = 3;

	Cij_l[0] = 2.0f;
	Cij_l[1] = -5.0f;
	Cij_l[2] = 7.0f;
	Cij_l[3] = 4.0f;
	Cij_l[4] = -1.0f;
	Cij_l[5] = 6.0f;
	Cij_l[6] = 3.0f;

	phaseI_setup(N, L, B_n, Ni_l, Nj_l);

	free(B_n);

	free(Basis_l);
	free(Ni_l);
	free(Nj_l);
	free(Cij_l);
	free(Xij_l);
}

void phaseI_Problem_3()
{
	unsigned int N = 8;
	float* B_n = malloc(N * sizeof(float));

	B_n[0] = 10.0f;
	B_n[1] = 20.0f;
	B_n[2] = 0.0f;
	B_n[3] = -5.0f;
	B_n[4] = 0.0f;
	B_n[5] = 0.0f;
	B_n[6] = -15.0f;
	B_n[7] = -10.0f;

	//put the check for B_n sums here

	unsigned L = 11;
	MCFBool* Basis_l = malloc(L * sizeof(MCFBool));
	unsigned int* Ni_l = malloc(L * sizeof(unsigned int));
	unsigned int* Nj_l = malloc(L * sizeof(unsigned int));
	float* Cij_l = malloc(L * sizeof(float));
	float* Xij_l = malloc(L * sizeof(float));

	Ni_l[0] = 0;
	Ni_l[1] = 1;
	Ni_l[2] = 1;
	Ni_l[3] = 1;
	Ni_l[4] = 2;
	Ni_l[5] = 2;
	Ni_l[6] = 3;
	Ni_l[7] = 4;
	Ni_l[8] = 4;
	Ni_l[9] = 5;
	Ni_l[10] = 6;

	Nj_l[0] = 3;
	Nj_l[1] = 0;
	Nj_l[2] = 2;
	Nj_l[3] = 5;
	Nj_l[4] = 3;
	Nj_l[5] = 4;
	Nj_l[6] = 6;
	Nj_l[7] = 5;
	Nj_l[8] = 6;
	Nj_l[9] = 7;
	Nj_l[10] = 7;

	Cij_l[0] = 2.0f;
	Cij_l[1] = 1.0f;
	Cij_l[2] = 0.0f;
	Cij_l[3] = 6.0f;
	Cij_l[4] = 1.0f;
	Cij_l[5] = 4.0f;
	Cij_l[6] = 5.0f;
	Cij_l[7] = 2.0f;
	Cij_l[8] = 7.0f;
	Cij_l[9] = 8.0f;
	Cij_l[10] = 9.0f;

	phaseI_setup(N, L, B_n, Ni_l, Nj_l);

	free(B_n);

	free(Basis_l);
	free(Ni_l);
	free(Nj_l);
	free(Cij_l);
	free(Xij_l);
}

void phaseI_Problem()
{
	unsigned int N = 4;
	float* B_n = malloc(N * sizeof(float));

	B_n[0] = 3.0f;
	B_n[1] = 3.0f;
	B_n[2] = 1.0f;
	B_n[3] = -7.0f;

	//put the check for B_n sums here

	unsigned L = 5;
	MCFBool* Basis_l = malloc(L * sizeof(MCFBool));
	unsigned int* Ni_l = malloc(L * sizeof(unsigned int));
	unsigned int* Nj_l = malloc(L * sizeof(unsigned int));
	float* Cij_l = malloc(L * sizeof(float));
	float* Xij_l = malloc(L * sizeof(float));

	Ni_l[0] = 0;
	Ni_l[1] = 0;
	Ni_l[2] = 1;
	Ni_l[3] = 2;
	Ni_l[4] = 2;

	Nj_l[0] = 1;
	Nj_l[1] = 2;
	Nj_l[2] = 3;
	Nj_l[3] = 1;
	Nj_l[4] = 3;


	Cij_l[0] = 3.0f;
	Cij_l[1] = 4.0f;
	Cij_l[2] = 2.0f;
	Cij_l[3] = -1.0f;
	Cij_l[4] = 6.0f;

	phaseI_setup(N, L, B_n, Ni_l, Nj_l);

	free(B_n);

	free(Basis_l);
	free(Ni_l);
	free(Nj_l);
	free(Cij_l);
	free(Xij_l);
}

int main(int argc, const char** argv)
{
	//phaseI_Problem_Example();  //checks out
	//phaseII_problem_example();  //check out
	//phaseI_Problem();  //something we wrong, get the original code that I wrote
	//phaseII_problem();  //checks out
	//phaseI_Problem_3();  //just click enter once.  something has gone wrong here
	//phaseII_problem_3();
	//WHERE IS phaseII_problem_4()?!?!?!?!?!?
	phaseII_problem_4();  //HANGS UP SHIT!!!

	printf("\nPress something to exit.\n");
	getchar();

	return 0;
}