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

typedef struct{
	Cord *desc_;
	int i_, j_;
} ArgCampoEff;




void inicia(), *calccampo(void* k), *calccampoeff(void* k), *calcvect(void* arg), Checa(int qnt, char*argsv[]), reseta(), *att(void* arg);
void escreve(int *i);
void Reset();
double random();
ArgCross *SetParam(Cord* vetor1, Cord* vetor2, Cord* vetorRes);

Cord ***s0, ***s, ***H, ***dsdt, ***pos, ***Heff;
Cord ***aux, *desc_, *r1, *r2, *r3, *r4;
Cord ***jj;
Cord *saida;
ArgCampoEff *agrseff;
ArgCampo *campoext;
ArgGrid *variavel_integra;
FILE *fp, *ani, *dados, *field, *posf;
Cord *vects[numthreads * 2];
ArgCampoEff *args;
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
	posf = fopen("pos.dat", "w+");
	s0 = malloc(sizeof(Cord) * nspinx * nspiny);
	aux = malloc(sizeof(Cord) * nspinx * nspiny);
	pos = malloc(sizeof(Cord) * nspinx * nspiny);
	s = malloc(sizeof(Cord) * nspinx * nspiny);
	jj = malloc(sizeof(Cord) * nspinx * nspiny);
	H = malloc(sizeof(Cord) * nspinx * nspiny);
	Heff = malloc(sizeof(Cord) * nspinx * nspiny);
	dsdt = malloc(sizeof(Cord) * nspinx * nspiny);
	r1 = malloc(sizeof(Cord));
	r2 = malloc(sizeof(Cord));
	r3 = malloc(sizeof(Cord));
	r4 = malloc(sizeof(Cord));
	saida = malloc(sizeof(Cord));
	campoext = malloc(sizeof(ArgCampo));
	variavel_integra = malloc(sizeof(ArgGrid));
	desc_ = malloc(sizeof(Cord));
	desc_->x = desc_->y = desc_->z = 0.0;
	args = malloc(sizeof(ArgCampoEff));
	ArgCampo init;
	init.t_ = 0.0;
	for(i = 0; i < nspinx; i++){
		s0[i] = malloc(sizeof(Cord) * nspiny);
		aux[i] = malloc(sizeof(Cord) * nspiny);
		pos[i] = malloc(sizeof(Cord) * nspiny);
		s[i] = malloc(sizeof(Cord) * nspiny);
		jj[i] = malloc(sizeof(Cord) * nspiny);
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
			s0[i][j]->x = 2.0 * random() - 1.0;
			s0[i][j]->y = 2.0 * random() - 1.0;
			s0[i][j]->z = 2.0 * random() - 1.0;
			double norm = sqrt(s0[i][j]->x * s0[i][j]->x + s0[i][j]->y * s0[i][j]->y + s0[i][j]->z * s0[i][j]->z);
			s0[i][j]->x = s0[i][j]->x / norm;
			s0[i][j]->y = s0[i][j]->y / norm;
			s0[i][j]->z = s0[i][j]->z / norm;
			/*s0[i][j]->x = 1.0 / sqrt(3.0);
			s0[i][j]->y = 1.0 / sqrt(3.0);
			s0[i][j]->z = 1.0 / sqrt(3.0);*/
			pos[i][j]->x = i;
			pos[i][j]->y = j;
			pos[i][j]->z = 0;
			jj[i][j]->x = 0.0;
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
	H[k_->i_][k_->j_]->z = 0.0;
	return NULL;
}

double random(){
	return (double)rand() / (double)RAND_MAX;
}

void* calccampoeff(void *k){
	ArgCampoEff *k_ = (ArgCampoEff*) k;
	Cord *hmmm;
	aux[k_->i_][k_->j_]->x += k_->desc_->x;
	aux[k_->i_][k_->j_]->y += k_->desc_->y;
	aux[k_->i_][k_->j_]->z += k_->desc_->z;
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
}


void* Calcdsdt(void* arg){
	ArgGrid *k_ = (ArgGrid*) arg;
	args->i_ = k_->i_;
	args->j_ = k_->j_;
	args->desc_ = desc_;
	calccampoeff((void*) args);
	calcvect((void*)SetParam(aux[k_->i_][k_->j_], jj[k_->i_][k_->j_], vects[1]));
	saida->x = Heff[k_->i_][k_->j_]->x + vects[1]->x;
	saida->y = Heff[k_->i_][k_->j_]->y + vects[1]->y;
	saida->z = Heff[k_->i_][k_->j_]->z + vects[1]->z;
	calcvect((void*)SetParam(aux[k_->i_][k_->j_], saida, vects[1]));
	r1->x = vects[1]->x * h / 2.0;
	r1->y = vects[1]->y * h / 2.0;
	r1->z = vects[1]->z * h / 2.0;
	reseta();

	args->desc_ = r1;
	calccampoeff((void*) args);
	calcvect((void*)SetParam(aux[k_->i_][k_->j_], jj[k_->i_][k_->j_], vects[2]));
	saida->x = Heff[k_->i_][k_->j_]->x + vects[2]->x;
	saida->y = Heff[k_->i_][k_->j_]->y + vects[2]->y;
	saida->z = Heff[k_->i_][k_->j_]->z + vects[2]->z;
	calcvect((void*)SetParam(aux[k_->i_][k_->j_], saida, vects[2]));
	r2->x = vects[2]->x * h / 2.0;
	r2->y = vects[2]->y * h / 2.0;
	r2->z = vects[2]->z * h / 2.0;
	reseta();

	args->desc_ = r2;
	calccampoeff((void*) args);
	calcvect((void*)SetParam(aux[k_->i_][k_->j_], jj[k_->i_][k_->j_], vects[3]));
	saida->x = Heff[k_->i_][k_->j_]->x + vects[3]->x;
	saida->y = Heff[k_->i_][k_->j_]->y + vects[3]->y;
	saida->z = Heff[k_->i_][k_->j_]->z + vects[3]->z;
	calcvect((void*)SetParam(aux[k_->i_][k_->j_], saida, vects[3]));
	r3->x = vects[3]->x * h;
	r3->y = vects[3]->y * h;
	r3->z = vects[3]->z * h;
	reseta();

	args->desc_ = r3;
	calccampoeff((void*) args);
	calcvect((void*)SetParam(aux[k_->i_][k_->j_], jj[k_->i_][k_->j_], vects[4]));
	saida->x = Heff[k_->i_][k_->j_]->x + vects[4]->x;
	saida->y = Heff[k_->i_][k_->j_]->y + vects[4]->y;
	saida->z = Heff[k_->i_][k_->j_]->z + vects[4]->z;
	calcvect((void*)SetParam(aux[k_->i_][k_->j_], saida, vects[4]));
	r4->x = vects[4]->x * h;
	r4->y = vects[4]->y * h;
	r4->z = vects[4]->z * h;
	reseta();
	dsdt[k_->i_][k_->j_]->x = (2.0 * r1->x + 4.0 * r2->x + 2.0 * r3->x + r4->x) / 6.0;
	dsdt[k_->i_][k_->j_]->y = (2.0 * r1->y + 4.0 * r2->y + 2.0 * r3->y + r4->y) / 6.0;
	dsdt[k_->i_][k_->j_]->z = (2.0 * r1->z + 4.0 * r2->z + 2.0 * r3->z + r4->z) / 6.0;
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
	s[k_->i_][k_->j_]->x = s0[k_->i_][k_->j_]->x + dsdt[k_->i_][k_->j_]->x;
	s[k_->i_][k_->j_]->y = s0[k_->i_][k_->j_]->y + dsdt[k_->i_][k_->j_]->y;
	s[k_->i_][k_->j_]->z = s0[k_->i_][k_->j_]->z + dsdt[k_->i_][k_->j_]->z;
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
		escreve(&t);
		Reset();
	}
}
int aux___ = 0;
void escreve(int *i){
	if(*i % corte == 0){
		for(int i_ = 0; i_ < nspinx; i_++){
			for(int j_ = 0; j_ < nspiny; j_++){
				fprintf(fp, "%.12f \t %.12f \t %.12f \t ", s[i_][j_]->x, s[i_][j_]->y, s[i_][j_]->z);
				fprintf(posf, "%.12f \t %.12f \t %.12f \t ", pos[i_][j_]->x, pos[i_][j_]->y, pos[i_][j_]->z);
			}
		}
		fprintf(fp, "\n");
		fprintf(posf, "\n");
		aux___+=1;
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
	printf("Starting system:\nn = %d\nlambda = %.5f\n(tmin, tmax) = (%.5f, %.5f)\nh = %.5f\ncut = %d\nRunge-Kutta 4\n", n, lambda, tmin, tmax, h, corte);
	fprintf(dados, "Starting system:\nn = %d\nlambda = %.5f\n(tmin, tmax) = (%.5f, %.5f)\nh = %.5f\ncut = %d\nRunge-Kutta 4\n", n, lambda, tmin, tmax, h, corte);
}