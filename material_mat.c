void calc_De_Iso_PSS_Mat(double *De, const MatData *tab)
{
	IsoMatData *p = (IsoMatData *)tab;
	double f = p->young / (1 - p->poisson * p->poisson);
	double *D[] = { De, De+3, De+6, };
	p->shar_mod=p->young / (1 + p->poisson) / 2;
	D[0][0] = D[1][1] = f;
	D[0][1] = D[1][0] = f * p->poisson;
	D[0][2] = D[1][2] = D[2][0] = D[2][1]= 0;
	D[2][2] = p->shar_mod;
}

void calc_De_Iso_PSN_Mat(double *De, const MatData *tab)
{
	IsoMatData *p = (IsoMatData *)tab;
	double f = p->young / (1 - 2 * p->poisson) / (1 + p->poisson);
	double *D[] = { De, De+3, De+6, };
	p->shar_mod=p->young / (1 + p->poisson) / 2;
	D[0][0] = D[1][1] = f * (1 - p->poisson);
	D[0][1] = D[1][0] = f * p->poisson;
	D[0][2] = D[1][2] = D[2][0] = D[2][1]= 0;
	D[2][2] = p->shar_mod;
}

void calc_De_Iso_AXSol_Mat(double *De, const MatData *tab)
{
	IsoMatData *p = (IsoMatData *)tab;
	double f = p->young / (1 - 2 * p->poisson) / (1 + p->poisson);
	double *D[] = { De, De+4, De+8, De+12, };
	int i;
	for (i = 0; i < 16; i++) De [i] = 0;
	p->shar_mod=p->young / (1 + p->poisson) / 2;
	D[0][0] = D[1][1] = D[2][2] = f * (1 - p->poisson);
	D[0][1] = D[0][2] = D[1][0] = D[1][2] = D[2][0] = D[2][1] = f * p->poisson;
	D[3][3] = p->shar_mod;
}

void calc_De_Iso_Solid_Mat(double *De, const MatData *tab)
{
	IsoMatData *p = (IsoMatData *)tab;
	double f = p->young / (1 - 2 * p->poisson) / (1 + p->poisson);
	double *D[] = { De, De+6, De+12, De+18, De+24, De+30, };
	int i;
	for (i = 0; i < 36; i++) De [i] = 0;
	p->shar_mod=p->young / (1 + p->poisson) / 2;
	D[0][0] = D[1][1] = D[2][2] = f * (1 - p->poisson);
	D[0][1] = D[0][2] = D[1][0] = D[1][2] = D[2][0] = D[2][1] = f * p->poisson;
	D[3][3] = D[4][4] = D[5][5] = p->shar_mod;
}
