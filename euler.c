#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#define numthreads 10

typedef struct{
	double x, y, z;
} Cord;

typedef struct{
	Cord *vetor1, *vetor2, *vetorRes;
} ArgCross;

typedef struct{
	double t_;
	int i_, j_;
} ArgCampo;

typedef struct{
	int i_, j_;
} ArgGrid;




void inicia(), *calccampo(void* k), *calccampoeff(void* k), *calcvect(void* arg), Checa(int qnt, char*argsv[]), reseta(), *att(void* arg);
void escreve(int *i);
void Reset();
ArgCross *SetParam(Cord* vetor1, Cord* vetor2, Cord* vetorRes);

Cord ***s0, ***s, ***H, ***dsdt, ***pos, ***Heff;
Cord ***jj;
Cord ***aux;
Cord *saida;
ArgCampo *campoext;
ArgGrid *variavel_integra;
FILE *fp, *ani, *dados, *field, *posf;
Cord *vects[numthreads * 2];
pthread_t threads[numthreads];
int i, j;
int n = 1000;
int nspinx = 30, nspiny = 30;
double tmax = 10.0;
double tmin = 0.0;
double lambda = 0.1;
double t, *taux = NULL, h;
double A = 1.0, B = 1.0, C = 1.0;
int corte = 10;


void inicia(){
	h = (tmax - tmin) / n;
	fp = fopen("out.dat", "w+");
	field = fopen("field.dat", "w+");
	posf = fopen("pos.dat", "w+");
	s0 = malloc(sizeof(Cord) * nspinx * nspiny);
	aux = malloc(sizeof(Cord) * nspinx * nspiny);
	pos = malloc(sizeof(Cord) * nspinx * nspiny);
	s = malloc(sizeof(Cord) * nspinx * nspiny);
	jj = malloc(sizeof(Cord) * nspinx * nspiny);
	H = malloc(sizeof(Cord) * nspinx * nspiny);
	Heff = malloc(sizeof(Cord) * nspinx * nspiny);
	dsdt = malloc(sizeof(Cord) * nspinx * nspiny);
	saida = malloc(sizeof(Cord));
	campoext = malloc(sizeof(ArgCampo));
	variavel_integra = malloc(sizeof(ArgGrid));
	ArgCampo init;
	init.t_ = 0.0;
	for(i = 0; i < nspinx; i++){
		s0[i] = malloc(sizeof(Cord) * nspiny);
		jj[i] = malloc(sizeof(Cord) * nspiny);
		aux[i] = malloc(sizeof(Cord) * nspiny);
		pos[i] = malloc(sizeof(Cord) * nspiny);
		s[i] = malloc(sizeof(Cord) * nspiny);
		H[i] = malloc(sizeof(Cord) * nspiny);
		Heff[i] = malloc(sizeof(Cord) * nspiny);
		dsdt[i] = malloc(sizeof(Cord) * nspiny);
		for(j = 0; j < nspiny; j++){
			init.i_ = i;
			init.j_ = j;
			s0[i][j] = malloc(sizeof(Cord));
			aux[i][j] = malloc(sizeof(Cord));
			pos[i][j] = malloc(sizeof(Cord));
			s[i][j] = malloc(sizeof(Cord));
			jj[i][j] = malloc(sizeof(Cord));
			H[i][j] = malloc(sizeof(Cord));
			dsdt[i][j] = malloc(sizeof(Cord));
			Heff[i][j] = malloc(sizeof(Cord));
			s0[i][j]->x = (double) (1.0) / sqrt(3.0);
			s0[i][j]->y = (double) (1.0) / sqrt(3.0);
			s0[i][j]->z = (double) (1.0) / sqrt(3.0);
			pos[i][j]->x = i;
			pos[i][j]->y = j;
			pos[i][j]->z = 0;
			jj[i][j]->x = 1.0;
			jj[i][j]->y = 0.0;
			jj[i][j]->z = 0.0;
			calccampo((void*)&init);
		}
	}
	for(i = 0; i < numthreads; i++){
		vects[i] = malloc(sizeof(Cord));
	}
	reseta();
}


void* calccampo(void *k){
	ArgCampo* k_ = (ArgCampo*) k;
	H[k_->i_][k_->j_]->x = 0.0;
	H[k_->i_][k_->j_]->y = 0.0;
	H[k_->i_][k_->j_]->z = 1.0;
	return NULL;
}

void* calccampoeff(void *k){
	ArgGrid *k_ = (ArgGrid*) k;
	Cord *hmmm;
	hmmm = malloc(sizeof(Cord));
	hmmm->x = hmmm->y = hmmm->z = 0;
	for(int i = -1; i <= 1; i++){
		int entx = k_->i_ + i;
		if(entx < 0){
			entx = nspinx - 1;
		}else if(entx > nspinx - 1){
			entx = 0;
		}
		for(int j = -1; j<= 1; j++){
			int enty = k_->j_ + j;
			if(enty < 0){
				enty = nspiny - 1;
			}else if(enty > nspiny - 1){
				enty = 0;
			}
			hmmm->x += aux[entx][enty]->x;
			hmmm->y += aux[entx][enty]->y;
			hmmm->z += aux[entx][enty]->z;
		}
	}
	hmmm->x -= aux[k_->i_][k_->j_]->x;
	hmmm->y -= aux[k_->i_][k_->j_]->y;
	hmmm->z -= aux[k_->i_][k_->j_]->z;
	hmmm->x += A * aux[k_->i_][k_->j_]->x + H[k_->i_][k_->j_]->x;
	hmmm->y += B * aux[k_->i_][k_->j_]->y + H[k_->i_][k_->j_]->y;
	hmmm->z += C * aux[k_->i_][k_->j_]->z + H[k_->i_][k_->j_]->z;
	calcvect((void*)SetParam(aux[k_->i_][k_->j_], hmmm, vects[0]));
	Heff[k_->i_][k_->j_]->x = hmmm->x - lambda * vects[0]->x; //inverte o sinal pq inverte no produto cross
	Heff[k_->i_][k_->j_]->y = hmmm->y - lambda * vects[0]->y;
	Heff[k_->i_][k_->j_]->z = hmmm->z - lambda * vects[0]->z;
	free(hmmm);
	return NULL;
}


void* Calcdsdt(void* arg){
	ArgGrid *k_ = (ArgGrid*) arg;
	calccampoeff(arg);
	calcvect((void*)SetParam(s0[k_->i_][k_->j_], jj[k_->i_][k_->j_], vects[1]));
	saida->x = Heff[k_->i_][k_->j_]->x + vects[1]->x;
	saida->y = Heff[k_->i_][k_->j_]->y + vects[1]->y;
	saida->z = Heff[k_->i_][k_->j_]->z + vects[1]->z;
	calcvect((void*)SetParam(s0[k_->i_][k_->j_], saida, vects[1]));
	dsdt[k_->i_][k_->j_]->x = vects[1]->x;
	dsdt[k_->i_][k_->j_]->y = vects[1]->y;
	dsdt[k_->i_][k_->j_]->z = vects[1]->z;
	reseta();
	return NULL;
}

void* calcvect(void *arg){
	ArgCross* k = (ArgCross*) arg;
	k->vetorRes->x = k->vetor1->y * k->vetor2->z - k->vetor1->z * k->vetor2->y;
	k->vetorRes->y = -(k->vetor1->x * k->vetor2->z - k->vetor1->z * k->vetor2->x);
	k->vetorRes->z = k->vetor1->x * k->vetor2->y - k->vetor1->y * k->vetor2->x;
	free(k);
}

ArgCross *SetParam(Cord* vetor1, Cord* vetor2, Cord* vetorRes){
	ArgCross* saida = malloc(sizeof(ArgCross));
	saida->vetor1 = vetor1;
	saida->vetor2 = vetor2;
	saida->vetorRes = vetorRes;
	return saida;
}

void reseta(){
	for(int i = 0; i < nspinx; i++){
		for(int j = 0; j < nspiny; j++){
			aux[i][j]->x = s0[i][j]->x;
			aux[i][j]->y = s0[i][j]->y;
			aux[i][j]->z = s0[i][j]->z;
		}
	}
}

void* att(void* arg){
	ArgGrid* k_ = (ArgGrid*) arg;
	s[k_->i_][k_->j_]->x = s0[k_->i_][k_->j_]->x + dsdt[k_->i_][k_->j_]->x * h;
	s[k_->i_][k_->j_]->y = s0[k_->i_][k_->j_]->y + dsdt[k_->i_][k_->j_]->y * h;
	s[k_->i_][k_->j_]->z = s0[k_->i_][k_->j_]->z + dsdt[k_->i_][k_->j_]->z * h;
}

void Reset(){
	for(int i = 0; i < nspinx; i++){
		for(int j = 0; j < nspiny; j++){
			s0[i][j]->x = s[i][j]->x;
			s0[i][j]->y = s[i][j]->y;
			s0[i][j]->z = s[i][j]->z;
		}
	}
}

void Integra(){
	for(int t = 0; t <= n; t++){
		campoext->t_ = t * h;
		for(int i = 0; i < nspinx; i++){
			for(int j = 0; j < nspiny; j++){
				campoext->i_ = i;
				campoext->j_ = j;
				variavel_integra->i_ = i;
				variavel_integra->j_ = j;
				calccampo((void*) campoext);
				Calcdsdt((void*)variavel_integra);
				att((void*)variavel_integra);
			}
		}
		Reset();
		escreve(&t);
	}
}

void escreve(int *i){
	if(*i % corte == 0){
		for(int i = 0; i < nspinx; i++){
			for(int j = 0; j < nspiny; j++){
				fprintf(fp, "%.12f \t %.12f \t %.12f \t ", s[i][j]->x, s[i][j]->y, s[i][j]->z);
				fprintf(field, "%.12f \t %.12f \t %.12f \t ", H[i][j]->x, H[i][j]->y, H[i][j]->z);
				fprintf(posf, "%.12f \t %.12f \t %.12f \t ", pos[i][j]->x, pos[i][j]->y, pos[i][j]->z);
			}
		}
		fprintf(fp, "\n");
		fprintf(field, "\n");
	}
}

void Checa(int qnt, char*argsv[]){
	if (qnt == 6){
		n = atoi(argsv[1]);
		tmin = atof(argsv[2]);
		tmax = atof(argsv[3]);
		lambda = atof(argsv[4]);
		corte = atoi(argsv[5]);
	}else if(qnt == 5){
		n = atoi(argsv[1]);
		tmin = atof(argsv[2]);
		tmax = atof(argsv[3]);
		lambda = atof(argsv[4]);
	}else if(qnt == 4){
		n = atoi(argsv[1]);
		tmin = atof(argsv[2]);
		tmax = atof(argsv[3]);
	}else if(qnt == 3){
		n = atoi(argsv[1]);
		tmax = atof(argsv[2]);
	}else if(qnt == 2){
		n = atoi(argsv[1]);
	}
}

void Dados(){
	dados = fopen("data.dat", "w+");
	printf("Starting system:\nn = %d\nlambda = %.5f\n(tmin, tmax) = (%.5f, %.5f)\nh = %.5f\ncut = %d\nEuler\n", n, lambda, tmin, tmax, h, corte);
	fprintf(dados, "Starting system:\nn = %d\nlambda = %.5f\n(tmin, tmax) = (%.5f, %.5f)\nh = %.5f\ncut = %d\nEuler\n", n, lambda, tmin, tmax, h, corte);
}