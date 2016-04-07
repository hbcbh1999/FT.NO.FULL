/*******************************************************************
 * 		CARTESIAN.c
 *******************************************************************/
/*
      To-Do: Line 392, the values of ilower and iupper 
      
*/

#include "solver.h"
#include "melting.h"


//-------------------------------------------------------------------
//		INC_STATE
//-------------------------------------------------------------------

void INC_STATE::setZero(void)
{
	T = 0.0;
	D = 0.0;
}


//----------------------------------------------------------------
//		RECTANGLE
//----------------------------------------------------------------

//RECTANGLE::RECTANGLE()
RECTANGLE::RECTANGLE(): index(-1), comp(-1)
{
}

void RECTANGLE::setCoords(
	double *crds,
	int dim)
{
	int i;
	for (i = 0; i < dim; ++i)
	    coords[i] = crds[i];
}
//--------------------------------------------------------------------------
// 		CARTESIAN
//--------------------------------------------------------------------------

CARTESIAN::~CARTESIAN()
{
}

//---------------------------------------------------------------
//	initMesh
// include the following parts
// 1) setup edges:		
// 2) setup cell_center
//---------------------------------------------------------------
void CARTESIAN::initMesh(void)
{
	int i,j,k, index;
	double crds[MAXD];
	int icoords[MAXD];
	int num_cells;
	int cell_index;
	RECTANGLE       rectangle;

	// init vertices,edges & cell_center
	if (debugging("trace")) printf("Entering initMesh()\n");
	FT_MakeGridIntfc(front);
	setDomain();

	num_cells = 1;
	for (i = 0; i < dim; ++i)
	{
	    num_cells *= (top_gmax[i] + 1);
	}
	cell_center.insert(cell_center.end(),num_cells,rectangle);
	
	// setup vertices
	// left to right, down to up
	switch (dim)
	{
	case 1:
	    for (j = 0; j <= top_gmax[1]; j++)
	    {
	    	crds[0] = top_L[0] + top_h[0]*i;
		index = d_index1d(i,top_gmax);
	    	cell_center[index].setCoords(crds,dim);
	    	cell_center[index].icoords[0] = i;
	    }
	case 2:
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	crds[0] = top_L[0] + top_h[0]*i;
	    	crds[1] = top_L[1] + top_h[1]*j;
		index = d_index2d(i,j,top_gmax);
	    	cell_center[index].setCoords(crds,dim);
	    	cell_center[index].icoords[0] = i;
	    	cell_center[index].icoords[1] = j;
	    }
	    break;
	case 3:
	    for (k = 0; k <= top_gmax[2]; k++)
	    for (j = 0; j <= top_gmax[1]; j++)
	    for (i = 0; i <= top_gmax[0]; i++)
	    {
	    	crds[0] = top_L[0] + top_h[0]*i;
	    	crds[1] = top_L[1] + top_h[1]*j;
	    	crds[2] = top_L[2] + top_h[2]*k;
		index = d_index3d(i,j,k,top_gmax);
	    	cell_center[index].setCoords(crds,dim);
	    	cell_center[index].icoords[0] = i;
	    	cell_center[index].icoords[1] = j;
	    	cell_center[index].icoords[2] = k;
	    }
	}
	setComponent();
	FT_FreeGridIntfc(front);
	if (debugging("trace")) printf("Leaving initMesh()\n");
}

void CARTESIAN::setComponent(void)
{
	int i;
	
	// cell center components
	for (i = 0; i < cell_center.size(); i++)
	{
	    cell_center[i].comp = 
	    		getComponent(cell_center[i].icoords);
	}
}

void CARTESIAN::setInitialCondition(void)
{
	int i;
	double coords[MAXD];
	INTERFACE *intfc = front->interf;
	POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	STATE *sl,*sr;

	FT_MakeGridIntfc(front);
        setDomain();

	/* Initialize states at the interface */
        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            if (negative_component(hs) == LIQUID_COMP)
                *sl = (STATE)eqn_params->T0[1];
            else if (negative_component(hs) == SOLID_COMP)
                *sl = (STATE)eqn_params->T0[0];
	    else
		*sl = 0.0;
            if (positive_component(hs) == LIQUID_COMP)
                *sr = (STATE)eqn_params->T0[1];
            else if (positive_component(hs) == SOLID_COMP)
                *sr = (STATE)eqn_params->T0[0];
	    else
		*sr = 0.0;
        }

	// cell_center
	for (i = 0; i < cell_center.size(); i++)
	{
	    getInitialState(cell_center[i]);
	    eqn_params->field->temperature[i] = cell_center[i].state.T;
	}
}	/* end setInitialCondition */

void CARTESIAN::setIndexMap(COMPONENT sub_comp)
{
	static boolean first = YES;
	int i,j,k,ic,index;
	int llbuf[MAXD],uubuf[MAXD];

	if (debugging("trace")) printf("Entering setIndexMap()\n");
	if (first)
	{
	    first = NO;
	    switch (dim)
	    {
	    case 1:
	    	FT_VectorMemoryAlloc((POINTER*)&i_to_I,top_gmax[0]+1,INT);
	    	break;
	    case 2:
	    	FT_MatrixMemoryAlloc((POINTER*)&ij_to_I,top_gmax[0]+1,top_gmax[1]+1,INT);
	    	break;
	    case 3:
	    	FT_TriArrayMemoryAlloc((POINTER*)&ijk_to_I,top_gmax[0]+1,top_gmax[1]+1,top_gmax[2]+1,
					INT);
	    	break;
	    }
	}

	index = 0;
	for (i = 0; i < dim; ++i)
	{
	    llbuf[i] = lbuf[i] != 0 ? lbuf[0] : 1;
	    uubuf[i] = ubuf[i] != 0 ? ubuf[0] : 1;
	}
	switch (dim)
	{
	case 1:
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index1d(i,top_gmax);
		if (cell_center[ic].comp == sub_comp)
		{
	    	    i_to_I[i] = index + ilower;
	    	    index++;
		}
		else
		    i_to_I[i] = -1;
	    }
	    FT_ParallelExchCellIndex(front,llbuf,uubuf,(POINTER)i_to_I);
	    break;
	case 2:
	    for (i = imin; i <= imax; i++)
	    for (j = jmin; j <= jmax; j++)
	    {
		ic = d_index2d(i,j,top_gmax);
		if (cell_center[ic].comp == sub_comp)
		{
	    	    ij_to_I[i][j] = index + ilower;
	    	    index++;
		}
		else
		    ij_to_I[i][j] = -1;
	    }
	    FT_ParallelExchCellIndex(front,llbuf,uubuf,(POINTER)ij_to_I);
	    break;
	case 3:
	    for (i = imin; i <= imax; i++)
	    for (j = jmin; j <= jmax; j++)
	    for (k = kmin; k <= kmax; k++)
	    {
		ic = d_index3d(i,j,k,top_gmax);
		if (cell_center[ic].comp == sub_comp)
		{
	    	    ijk_to_I[i][j][k] = index + ilower;
	    	    index++;
		}
		else
		    ijk_to_I[i][j][k] = -1;
	    }
	    FT_ParallelExchCellIndex(front,llbuf,uubuf,(POINTER)ijk_to_I);
	    break;
	}
	if (debugging("trace")) printf("Leaving setIndexMap()\n");

}	/* end setIndexMap */

void CARTESIAN::computeAdvection(COMPONENT sub_comp)
{
	if (debugging("trace")) printf("Entering computeAdvection()\n");

	switch (sub_comp)
	{
	case LIQUID_COMP:
	    eqn_params->D = eqn_params->k[1]/eqn_params->rho[1]
				/eqn_params->Cp[1];
	    break;
	case SOLID_COMP:
	    eqn_params->D = eqn_params->k[0]/eqn_params->rho[0]
				/eqn_params->Cp[0];
	    break;
	}

	switch (eqn_params->num_scheme)
	{
	case UNSPLIT_EXPLICIT:
	    return computeAdvectionExplicit(sub_comp);
	case UNSPLIT_IMPLICIT:
	    return computeAdvectionImplicit(sub_comp);
	case CRANK_NICOLSON:
	    return computeAdvectionCN(sub_comp);
	}
}
    
void CARTESIAN::computeAdvectionCN(COMPONENT sub_comp)
{
	int i,j,k,l,m,ic,icn,I,I_nb,icoords[MAXD];
	int gmin[MAXD],ipn[MAXD];
	double crx_coords[MAXD];
	double T0,T_nb,D,lambda,coeff,coeff_nb,rhs;
	COMPONENT comp;
	PETSc solver;
	double *x;
        int num_iter = 0;
        double rel_residual = 0;
        boolean fr_crx_grid_seg;
        const GRID_DIRECTION dir[3][2] =
                {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};

	start_clock("computeAdvectionCN");
	if (debugging("trace")) printf("Entering computeAdvectionCN()\n");

	for (i = 0; i < dim; ++i) gmin[i] = 0;

	D = eqn_params->D;
	setIndexMap(sub_comp);
	if (debugging("trace")) 
	{
	    int domain_size = 1;
	    printf("ilower = %d  iupper = %d\n",ilower,iupper);
	    for (i = 0; i < dim; ++i)
		domain_size *= (imax-imin+1);
	    printf("domain_size = %d\n",domain_size);
	}
        
	start_clock("set_coefficients");
	switch(dim)
	{
	case 1:
	    solver.Create(ilower, iupper-1, 3, 0);
	    for (i = imin; i <= imax; ++i)
	    {
	    	icoords[0] = i;
	    	ic = d_index1d(i,top_gmax);
	    	comp = top_comp[ic];
	    	if (comp != sub_comp) 
	    	    continue;
		I = i_to_I[i];
                T0 = cell_center[ic].state.T;
		rhs = T0;
		coeff = 1.0;
		for (l = 0; l < dim; ++l)
		{
		    lambda = D*m_dt/sqr(top_h[l]);
		    coeff += lambda;
                    for (m = 0; m < 2; ++m)
                    {
                        next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                        icn = d_index1d(ipn[0],top_gmax);
			I_nb = i_to_I[ipn[0]];
			coeff_nb = -0.5*lambda;
                        fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                icoords,dir[l][m],comp,temperature_of_state,
                                &T_nb,crx_coords);
                        if (!fr_crx_grid_seg)
                        {
			    solver.Add_A(I,I_nb,coeff_nb);
			    T_nb = cell_center[icn].state.T;
			    rhs -= coeff_nb*(T_nb - T0);
                        }
			else
			{
			    rhs -= coeff_nb*(2.0*T_nb - T0);
			}
                    }
		}
		solver.Add_A(I,I,coeff);
		solver.Add_b(I,rhs);
	    }
	    break;
	case 2:
	    solver.Create(ilower, iupper-1, 5, 0);
	    for (i = imin; i <= imax; ++i)
	    for (j = jmin; j <= jmax; ++j)
	    {
	    	icoords[0] = i;
	    	icoords[1] = j;
	    	ic = d_index2d(i,j,top_gmax);
	    	comp = top_comp[ic];
		I = ij_to_I[i][j];
	    	if (comp != sub_comp) 
	    	    continue;
                T0 = cell_center[ic].state.T;
		rhs = T0;
		coeff = 1.0;
		for (l = 0; l < dim; ++l)
		{
		    lambda = D*m_dt/sqr(top_h[l]);
		    coeff += lambda;
                    for (m = 0; m < 2; ++m)
                    {
                        next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                        icn = d_index2d(ipn[0],ipn[1],top_gmax);
			I_nb = ij_to_I[ipn[0]][ipn[1]];
			coeff_nb = -0.5*lambda;
                        fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                icoords,dir[l][m],comp,temperature_of_state,
                                &T_nb,crx_coords);
                        if (!fr_crx_grid_seg)
                        {
			    solver.Add_A(I,I_nb,coeff_nb);
			    T_nb = cell_center[icn].state.T;
			    rhs -= coeff_nb*(T_nb - T0);
                        }
			else
			{
			    rhs -= coeff_nb*(2.0*T_nb - T0);
			}
                    }
		}
		solver.Add_A(I,I,coeff);
		solver.Add_b(I,rhs);
	    }
	    break;
	case 3:
	    solver.Create(ilower, iupper-1, 7, 0);
	    for (i = imin; i <= imax; ++i)
	    for (j = jmin; j <= jmax; ++j)
	    for (k = kmin; k <= kmax; ++k)
	    {
	    	icoords[0] = i;
	    	icoords[1] = j;
	    	icoords[2] = k;
	    	ic = d_index3d(i,j,k,top_gmax);
	    	comp = top_comp[ic];
		I = ijk_to_I[i][j][k];
	    	if (comp != sub_comp) 
	    	    continue;
                T0 = cell_center[ic].state.T;
		rhs = T0;
		coeff = 1.0;
		for (l = 0; l < dim; ++l)
		{
		    lambda = D*m_dt/sqr(top_h[l]);
		    coeff += lambda;
                    for (m = 0; m < 2; ++m)
                    {
                        next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                        icn = d_index3d(ipn[0],ipn[1],ipn[2],top_gmax);
			I_nb = ijk_to_I[ipn[0]][ipn[1]][ipn[2]];
			coeff_nb = -0.5*lambda;
                        fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                icoords,dir[l][m],comp,temperature_of_state,
                                &T_nb,crx_coords);
                        if (!fr_crx_grid_seg)
                        {
			    solver.Add_A(I,I_nb,coeff_nb);
			    T_nb = cell_center[icn].state.T;
			    rhs -= coeff_nb*(T_nb - T0);
                        }
			else
			{
			    rhs -= coeff_nb*(2.0*T_nb - T0);
			}
                    }
		}
		solver.Add_A(I,I,coeff);
		solver.Add_b(I,rhs);
	    }
	    break;
	}
	stop_clock("set_coefficients");
	start_clock("petsc_solve");
	solver.SetMaxIter(500);   
	solver.SetTol(1e-8);   
	solver.Solve();

	if (debugging("PETSc"))
	{
	    (void) printf("CARTESIAN::computeAdvectionCN: "
	       		"num_iter = %d, rel_residual = %le \n", 
			num_iter, rel_residual);
	}

	FT_VectorMemoryAlloc((POINTER*)&x,cell_center.size(),sizeof(double));
	solver.Get_x(x);
	stop_clock("petsc_solve");

	start_clock("scatter_data");
	switch (dim)
	{
        case 1:
            for (i = imin; i <= imax; i++)
            {
	    	I = i_to_I[i];
	    	ic = d_index1d(i,top_gmax);
	    	comp = cell_center[ic].comp;
	    	if (comp == sub_comp)
	    	    array[ic] = x[I-ilower];
	    	else
	    	    array[ic] = 0.0;
	    }
	    break;
        case 2:
            for (i = imin; i <= imax; i++)
            for (j = jmin; j <= jmax; j++)
            {
	    	I = ij_to_I[i][j];
	    	ic = d_index2d(i,j,top_gmax);
	    	comp = cell_center[ic].comp;
	    	if (comp == sub_comp)
	    	    array[ic] = x[I-ilower];
	    	else
	    	    array[ic] = 0.0;
	    }
	    break;
        case 3:
            for (i = imin; i <= imax; i++)
            for (j = jmin; j <= jmax; j++)
            for (k = kmin; k <= kmax; k++)
            {
	    	I = ijk_to_I[i][j][k];
	    	ic = d_index3d(i,j,k,top_gmax);
	    	comp = cell_center[ic].comp;
	    	if (comp == sub_comp)
	    	    array[ic] = x[I-ilower];
	    	else
	    	    array[ic] = 0.0;
	    }
	    break;
	}
	scatMeshArray();
	switch (dim)
	{
        case 1:
            for (i = 0; i <= top_gmax[0]; ++i)
            {
		ic = d_index1d(i,top_gmax);
	    	comp = cell_center[ic].comp;
	    	if (comp == sub_comp)
		    cell_center[ic].state.T = array[ic];
	    }
	    break;
        case 2:
            for (i = 0; i <= top_gmax[0]; ++i)
            for (j = 0; j <= top_gmax[1]; ++j)
            {
		ic = d_index2d(i,j,top_gmax);
	    	comp = cell_center[ic].comp;
	    	if (comp == sub_comp)
		    cell_center[ic].state.T = array[ic];
	    }
	    break;
        case 3:
            for (i = 0; i <= top_gmax[0]; ++i)
            for (j = 0; j <= top_gmax[1]; ++j)
            for (k = 0; k <= top_gmax[2]; ++k)
            {
		ic = d_index3d(i,j,k,top_gmax);
	    	comp = cell_center[ic].comp;
	    	if (comp == sub_comp)
		    cell_center[ic].state.T = array[ic];
	    }
	    break;
	}
	stop_clock("scatter_data");

	if (debugging("trace")) printf("Leaving computeAdvectionCN()\n");
	stop_clock("computeAdvectionCN");
}	/* end computeAdvectionCN */

void CARTESIAN::computeAdvectionImplicit(COMPONENT sub_comp)
{
	int i,j,k,l,m,ic,icn,I,I_nb,icoords[MAXD];
	int gmin[MAXD],ipn[MAXD];
	double crx_coords[MAXD];
	double T0,T_nb,D,lambda,coeff,coeff_nb,rhs;
	COMPONENT comp;
	PETSc solver;
	double *x;
        int num_iter = 0;
        double rel_residual = 0;
        boolean fr_crx_grid_seg;
        const GRID_DIRECTION dir[3][2] =
                {{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};

	start_clock("computeAdvectionImplicit");
	if (debugging("trace")) printf("Entering computeAdvectionImplicit()\n");

	for (i = 0; i < dim; ++i) gmin[i] = 0;

	D = eqn_params->D;
	setIndexMap(sub_comp);
	if (debugging("trace")) 
	{
	    int domain_size = 1;
	    printf("ilower = %d  iupper = %d\n",ilower,iupper);
	    for (i = 0; i < dim; ++i)
		domain_size *= (imax-imin+1);
	    printf("domain_size = %d\n",domain_size);
	}
        
	start_clock("set_coefficients");
	switch(dim)
	{
	case 1:
	    solver.Create(ilower, iupper-1, 3, 0);
	    for (i = imin; i <= imax; ++i)
	    {
	    	icoords[0] = i;
	    	ic = d_index1d(i,top_gmax);
	    	comp = top_comp[ic];
	    	if (comp != sub_comp) 
	    	    continue;
		I = i_to_I[i];
                T0 = cell_center[ic].state.T;
		rhs = T0;
		coeff = 1.0;
		for (l = 0; l < dim; ++l)
		{
		    lambda = D*m_dt/sqr(top_h[l]);
		    coeff += 2.0*lambda;
                    for (m = 0; m < 2; ++m)
                    {
                        next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                        icn = d_index1d(ipn[0],top_gmax);
			I_nb = i_to_I[ipn[0]];
			coeff_nb = -lambda;
                        fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                icoords,dir[l][m],comp,temperature_of_state,
                                &T_nb,crx_coords);
                        if (!fr_crx_grid_seg)
                        {
			    solver.Add_A(I,I_nb,coeff_nb);
                        }
			else
			{
			    rhs -= coeff_nb*T_nb;
			}
                    }
		}
		solver.Add_A(I,I,coeff);
		solver.Add_b(I,rhs);
	    }
	    break;
	case 2:
	    solver.Create(ilower, iupper-1, 5, 0);
	    for (i = imin; i <= imax; ++i)
	    for (j = jmin; j <= jmax; ++j)
	    {
	    	icoords[0] = i;
	    	icoords[1] = j;
	    	ic = d_index2d(i,j,top_gmax);
	    	comp = top_comp[ic];
		I = ij_to_I[i][j];
	    	if (comp != sub_comp) 
	    	    continue;
                T0 = cell_center[ic].state.T;
		rhs = T0;
		coeff = 1.0;
		for (l = 0; l < dim; ++l)
		{
		    lambda = D*m_dt/sqr(top_h[l]);
		    coeff += 2.0*lambda;
                    for (m = 0; m < 2; ++m)
                    {
                        next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                        icn = d_index2d(ipn[0],ipn[1],top_gmax);
			I_nb = ij_to_I[ipn[0]][ipn[1]];
			coeff_nb = -lambda;
                        fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                icoords,dir[l][m],comp,temperature_of_state,
                                &T_nb,crx_coords);
                        if (!fr_crx_grid_seg)
                        {
			    solver.Add_A(I,I_nb,coeff_nb);
                        }
			else
			{
			    rhs -= coeff_nb*T_nb;
			}
                    }
		}
		solver.Add_A(I,I,coeff);
		solver.Add_b(I,rhs);
	    }
	    break;
	case 3:
	    solver.Create(ilower, iupper-1, 7, 0);
	    for (i = imin; i <= imax; ++i)
	    for (j = jmin; j <= jmax; ++j)
	    for (k = kmin; k <= kmax; ++k)
	    {
	    	icoords[0] = i;
	    	icoords[1] = j;
	    	icoords[2] = k;
	    	ic = d_index3d(i,j,k,top_gmax);
	    	comp = top_comp[ic];
		I = ijk_to_I[i][j][k];
	    	if (comp != sub_comp) 
	    	    continue;
                T0 = cell_center[ic].state.T;
		rhs = T0;
		coeff = 1.0;
		for (l = 0; l < dim; ++l)
		{
		    lambda = D*m_dt/sqr(top_h[l]);
		    coeff += 2.0*lambda;
                    for (m = 0; m < 2; ++m)
                    {
                        next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                        icn = d_index3d(ipn[0],ipn[1],ipn[2],top_gmax);
			I_nb = ijk_to_I[ipn[0]][ipn[1]][ipn[2]];
			coeff_nb = -lambda;
                        fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                icoords,dir[l][m],comp,temperature_of_state,
                                &T_nb,crx_coords);
                        if (!fr_crx_grid_seg)
                        {
			    solver.Add_A(I,I_nb,coeff_nb);
                        }
			else
			{
			    rhs -= coeff_nb*T_nb;
			}
                    }
		}
		solver.Add_A(I,I,coeff);
		solver.Add_b(I,rhs);
	    }
	    break;
	}
	stop_clock("set_coefficients");
	start_clock("petsc_solve");
	solver.SetMaxIter(500);   
	solver.SetTol(1e-6);   
	solver.Solve();

	if (debugging("PETSc"))
	{
	    (void) printf("CARTESIAN::computeAdvectionImplicit: "
	       		"num_iter = %d, rel_residual = %le \n", 
			num_iter, rel_residual);
	}

	FT_VectorMemoryAlloc((POINTER*)&x,cell_center.size(),sizeof(double));
	solver.Get_x(x);
	stop_clock("petsc_solve");

	start_clock("scatter_data");
	switch (dim)
	{
        case 1:
            for (i = imin; i <= imax; i++)
            {
	    	I = i_to_I[i];
	    	ic = d_index1d(i,top_gmax);
	    	comp = cell_center[ic].comp;
	    	if (comp == sub_comp)
	    	    array[ic] = x[I-ilower];
	    }
	    break;
        case 2:
            for (i = imin; i <= imax; i++)
            for (j = jmin; j <= jmax; j++)
            {
	    	I = ij_to_I[i][j];
	    	ic = d_index2d(i,j,top_gmax);
	    	comp = cell_center[ic].comp;
	    	if (comp == sub_comp)
	    	    array[ic] = x[I-ilower];
	    }
	    break;
        case 3:
            for (i = imin; i <= imax; i++)
            for (j = jmin; j <= jmax; j++)
            for (k = kmin; k <= kmax; k++)
            {
	    	I = ijk_to_I[i][j][k];
	    	ic = d_index3d(i,j,k,top_gmax);
	    	comp = cell_center[ic].comp;
	    	if (comp == sub_comp)
	    	    array[ic] = x[I-ilower];
	    }
	    break;
	}
	scatMeshArray();
	switch (dim)
	{
        case 1:
            for (i = 0; i <= top_gmax[0]; ++i)
            {
		ic = d_index1d(i,top_gmax);
	    	comp = cell_center[ic].comp;
	    	if (comp == sub_comp)
		    cell_center[ic].state.T = array[ic];
	    }
	    break;
        case 2:
            for (i = 0; i <= top_gmax[0]; ++i)
            for (j = 0; j <= top_gmax[1]; ++j)
            {
		ic = d_index2d(i,j,top_gmax);
	    	comp = cell_center[ic].comp;
	    	if (comp == sub_comp)
		    cell_center[ic].state.T = array[ic];
	    }
	    break;
        case 3:
            for (i = 0; i <= top_gmax[0]; ++i)
            for (j = 0; j <= top_gmax[1]; ++j)
            for (k = 0; k <= top_gmax[2]; ++k)
            {
		ic = d_index3d(i,j,k,top_gmax);
	    	comp = cell_center[ic].comp;
	    	if (comp == sub_comp)
		    cell_center[ic].state.T = array[ic];
	    }
	    break;
	}
	stop_clock("scatter_data");

	if (debugging("trace")) printf("Leaving computeAdvectionImplicit()\n");
	stop_clock("computeAdvectionImplicit");
}	/* end computeAdvectionImplicit */


void CARTESIAN::computeSourceTerm(double *coords, double t, INC_STATE &state) 
{
	//computeSourceTerm(coords, state);
}

void CARTESIAN::getInitialState(RECTANGLE &cell) 
{
	COMPONENT c;

	c = getComponent(cell.icoords);
	if (c == LIQUID_COMP)
	    cell.state.T = eqn_params->T0[1];
	else if (c == SOLID_COMP)
	    cell.state.T = eqn_params->T0[0];
	else
	    cell.state.T = 0.0;
}

// for initial condition: 
// 		setInitialCondition();	
// this function should be called before solve()
// for the source term of the momentum equation: 	
// 		computeSourceTerm();
void CARTESIAN::solve(double dt)
{

	if (debugging("trace")) printf("Entering solve()\n");
	m_dt = dt;
	start_clock("solve");

	setDomain();
        if (debugging("trace")) printf("Passing setDomain()\n");

	setComponent();
	if (debugging("trace")) printf("Passing setComponent()\n");

	setGlobalIndex(LIQUID_COMP);
	if (debugging("trace")) printf("Passing liquid setGlobalIndex()\n");
	computeAdvection(LIQUID_COMP);
	if (debugging("trace")) printf("Passing liquid computeAdvection()\n");

	setGlobalIndex(SOLID_COMP);
	if (debugging("trace")) printf("Passing solid setGlobalIndex()\n");
	computeAdvection(SOLID_COMP);
	if (debugging("trace")) printf("Passing solid computeAdvection()\n");

	copyTemp();
	if (debugging("trace")) printf("Passing copyTemp()\n");
	
	setAdvectionDt();
	if (debugging("trace")) printf("Passing setAdvectionDt()\n");

	setBoundary();
	if (debugging("trace")) printf("Passing setBoundary()\n");

	copyMeshStates();
        if (debugging("trace")) printf("Passing copySolute()\n");

	stop_clock("solve");
	if (debugging("trace")) printf("Leaving solve()\n");
}


void CARTESIAN::setAdvectionDt()
{
	double D,Dl,Ds;

	eqn_params = (PARAMS*)front->extra2;
	Dl = eqn_params->k[1]/eqn_params->rho[1]/eqn_params->Cp[1];
	Ds = eqn_params->k[0]/eqn_params->rho[0]/eqn_params->Cp[0];
	D = std::max(Dl,Ds);

	if (eqn_params->num_scheme == UNSPLIT_EXPLICIT)
	{
	    m_dt = 0.5*sqr(hmin)/D/(double)dim;
	}
	else
	{
	    m_dt = 0.5*hmin/D/(double)dim;
	}
	if (debugging("trace"))
	{
	    printf("In setAdvectionDt: m_dt = %24.18g\n",m_dt);
	}
}	/* end setAdvectionDt */

void CARTESIAN::getVelocity(double *p, double *U)
{
	// locate the point
	int icoords[MAXD];
	int i,j,k;
	double c1[MAXD], c2[MAXD], denominator;

	if (!rect_in_which(p,icoords,top_grid))
	{
	    for (i=0; i<2; i++)
	    {
	    	U[i] = 0.0;
	    }
	    return;
	}

	switch (dim)
	{
	case 2:
	    /*
	    int index11, index21, index12, index22;
	    i = icoords[0];
	    j = icoords[1];
	    index11 = d_index2d(i,j,top_gmax);
	    index12 = d_index2d(i,j+1,top_gmax);
	    index21 = d_index2d(i+1,j,top_gmax);
	    index22 = d_index2d(i+1,j+1,top_gmax);

	    getRectangleCenter(index11, c1);
	    getRectangleCenter(index22, c2);

	    denominator = (c2[0]-c1[0])*(c2[1]-c1[1]);
	    for (i=0; i<2; i++)
	    {
		U[i] =  cell_center[index11].state.m_U[i]/denominator
				*(c2[0]-p[0])*(c2[1]-p[1]) 
			  + cell_center[index21].state.m_U[i]/denominator 
			  	*(p[0]-c1[0])*(c2[1]-p[1])
			  + cell_center[index12].state.m_U[i]/denominator
			  	*(c2[0]-p[0])*(p[1]-c1[1])
			  + cell_center[index22].state.m_U[i]/denominator
			  	*(p[0]-c1[0])*(p[1]-c1[1]);
	    }
	    */
	    break;
	}
}

void CARTESIAN::getRectangleIndex(int index, int &i, int &j)
{
	i = cell_center[index].icoords[0];
	j = cell_center[index].icoords[1];
}

void CARTESIAN::getRectangleIndex(int index, int &i, int &j, int &k)
{
	i = cell_center[index].icoords[0];
	j = cell_center[index].icoords[1];
	k = cell_center[index].icoords[2];
}


int CARTESIAN::getRectangleComponent(int index)
{	
	return getComponent(cell_center[index].icoords);
}

void CARTESIAN::getRectangleCenter(
	int index, 
	double *coords)
{
	int i;
	for (i = 0; i < dim; ++i)
	    coords[i] = cell_center[index].coords[i];
}

void CARTESIAN::getRectangleCenter(
	int index0, 
	int index1, 
	double *coords)
{
	int i;
	for (i = 0; i < dim; ++i)
	{
	    coords[i] = 0.5*(cell_center[index0].coords[i] +
	    		     cell_center[index1].coords[i]);
	}
}

int CARTESIAN::getComponent(
	double *p)
{
	return component(p,front->interf);
}

int CARTESIAN::getComponent(
	int *icoords)
{
	int index;
	switch (dim)
	{
	case 1:
	    index = d_index1d(icoords[0],top_gmax);
	    return top_comp[index];
	case 2:
	    index = d_index2d(icoords[0],icoords[1],top_gmax);
	    return top_comp[index];
	case 3:
	    index = d_index3d(icoords[0],icoords[1],icoords[2],top_gmax);
	    return top_comp[index];
	}
}

void CARTESIAN::save(char *filename)
{
	
	RECT_GRID *rect_grid = front->rect_grid;
	INTERFACE *intfc    = front->interf;
		
	int i, j;
	int xmax = rect_grid->gmax[0];
	int ymax = rect_grid->gmax[1];
	double x, y;
	
	FILE *hfile = fopen(filename, "w");
	if(hfile==NULL)
	{
		printf("\n can't open %s in "
		       "SaveAsTecplot_rect_grid_and_interface().", filename);
		exit(0);
	}
	
	// secondly print out the interface
		
	if(exists_interface(intfc))
	{
	    CURVE		**cur;
	    CURVE		*curve;
	    BOND		*bond;
			
	    for(cur=intfc->curves; cur && *cur; cur++)
	    {
		curve = *cur;
		fprintf(hfile, "ZONE I=%d J=%d F=POINT \n", 
				curve->num_points, 1);
		bond=curve->first;
		fprintf(hfile, "%.4f %.4f \n",bond->start->_coords[0], 
				bond->start->_coords[1]);
		for(bond=curve->first; bond!=NULL; bond=bond->next)
		    fprintf(hfile, "%.4f %.4f \n",bond->end->_coords[0], 
		    		bond->end->_coords[1]);
		}					
	}		
	fclose(hfile);
}

CARTESIAN::CARTESIAN(Front &front):front(&front)
{
}

void CARTESIAN::makeGridIntfc()
{
	static boolean first = YES;
	INTERFACE *grid_intfc;
	Table *T;
	static double *temperature;
	int i;

	FT_MakeGridIntfc(front);
	if (debugging("trace")) printf("Passed FT_MakeGridIntfc()\n");

	grid_intfc = front->grid_intfc;
	top_grid = &topological_grid(grid_intfc);
	lbuf = front->rect_grid->lbuf;
	ubuf = front->rect_grid->ubuf;
	top_gmax = top_grid->gmax;
	top_L = top_grid->L;
	top_U = top_grid->U;
	top_h = top_grid->h;
	dim = grid_intfc->dim;
	T = table_of_interface(grid_intfc);
	top_comp = T->components;
	eqn_params = (PARAMS*)front->extra2;
	hmin = top_h[0];
	for (i = 1; i < dim; ++i)
	    if (hmin > top_h[i]) hmin = top_h[i];

	switch (dim)
	{
	case 1:
            if (first)
            {
                FT_VectorMemoryAlloc((POINTER*)&array,top_gmax[0]+1,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&temperature,top_gmax[0]+1,FLOAT);
                first = NO;
            }	
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    eqn_params->field->temperature = temperature;
	    eqn_params = (PARAMS*)front->extra2;
	    break;
	case 2:
	    if (first)
	    {
	    	FT_VectorMemoryAlloc((POINTER*)&array,(top_gmax[0]+1)*(top_gmax[1]+1),FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&temperature,(top_gmax[0]+1)*(top_gmax[1]+1),FLOAT);
	    	first = NO;
	    }
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
	    eqn_params->field->temperature = temperature;
	    break;
	case 3:
	    if (first)
	    {
	    	FT_VectorMemoryAlloc((POINTER*)&array,
			(top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1),FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&temperature,
			(top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1),FLOAT);
	    	first = NO;
	    }
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
	    kmin = (lbuf[2] == 0) ? 1 : lbuf[2];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
	    kmax = (ubuf[2] == 0) ? top_gmax[2] - 1 : top_gmax[2] - ubuf[2];
	    eqn_params->field->temperature = temperature;
	    break;
	}
}	/* end makeGridIntfc */

void CARTESIAN::deleteGridIntfc()
{
	FT_FreeGridIntfc(front);
}

void CARTESIAN::scatMeshArray()
{
	FT_ParallelExchGridArrayBuffer(array,front);
}

void CARTESIAN::setGlobalIndex(COMPONENT sub_comp)
{
	int i,j,k,ic;
	int num_nodes = pp_numnodes();
	int myid = pp_mynode();
	static boolean first = YES;

	if (first)
	{
	    first = NO;
	    FT_VectorMemoryAlloc((POINTER*)&n_dist,num_nodes,sizeof(int));
	}
	NLblocks = 0;
	switch (dim)
	{
	case 1:
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index1d(i,top_gmax);
		if (cell_center[ic].comp != sub_comp) continue;
		NLblocks++;
	    }
	    break;
	case 2:
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index2d(i,j,top_gmax);
		if (cell_center[ic].comp != sub_comp) continue;
		NLblocks++;
	    }
	    break;
	case 3:
	    for (k = kmin; k <= kmax; k++)
	    for (j = jmin; j <= jmax; j++)
	    for (i = imin; i <= imax; i++)
	    {
		ic = d_index3d(i,j,k,top_gmax);
		if (cell_center[ic].comp != sub_comp) continue;
		NLblocks++;
	    }
	    break;
	}

	for (i = 0; i < num_nodes; ++i) n_dist[i] = 0;
	n_dist[myid] = NLblocks;
	pp_global_imax(n_dist,num_nodes);
	ilower = 0;
        iupper = n_dist[0];

        for (i = 1; i <= myid; i++)
        {
            ilower += n_dist[i-1];
            iupper += n_dist[i];
        }	
}


void CARTESIAN::copyTemp()
{
	int i,j,k,index;
	double *temperature = eqn_params->field->temperature;

	switch (dim)
	{
	case 1:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
	    	index  = d_index1d(i,top_gmax);
	    	temperature[index] = cell_center[index].state.T;
	    }
	    break;
	case 2:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    {
	    	index  = d_index2d(i,j,top_gmax);
	    	temperature[index] = cell_center[index].state.T;
	    }
	    break;
	case 3:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (k = 0; k <= top_gmax[2]; ++k)
	    {
	    	index  = d_index3d(i,j,k,top_gmax);
	    	temperature[index] = cell_center[index].state.T;
	    }
	    break;
	}
}	/* end copyTemp */

void CARTESIAN::checkStates()
{
}

void CARTESIAN::printFrontInteriorState(char *out_name)
{
	int i,j,k,l,index;
	char filename[100];
	FILE *outfile;
	INTERFACE *intfc = front->interf;
	STATE *sl,*sr;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;

	sprintf(filename,"%s/state.ts%s",out_name,right_flush(front->step,7));
#if defined(__MPI__)
        sprintf(filename,"%s-nd%s",filename,right_flush(pp_mynode(),4));
#endif /* defined(__MPI__) */
	outfile = fopen(filename,"w");
	
        /* Initialize states at the interface */
        fprintf(outfile,"Interface states:\n");
        int count = 0;
        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            count++;
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            fprintf(outfile,"%24.18g %24.18g\n",(double)*sl,(double)*sr);
        }

	fprintf(outfile,"\nInterior states:\n");
	switch (dim)
	{
	case 1:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
	        fprintf(outfile,"%24.18g\n",cell_center[index].state.T);
	    }
	    break;
	case 2:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    {
		index = d_index2d(i,j,top_gmax);
	        fprintf(outfile,"%24.18g\n",cell_center[index].state.T);
	    }
	    break;
	case 3:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (k = 0; k <= top_gmax[2]; ++k)
	    {
		index = d_index3d(i,j,k,top_gmax);
	        fprintf(outfile,"%24.18g\n",
			cell_center[index].state.T);
	    }
	    break;
	}
	fclose(outfile);
}

void CARTESIAN::readFrontInteriorState(char *restart_name)
{
	FILE *infile;
	int i,j,k,l,index;
	INTERFACE *intfc = front->interf;
	STATE *sl,*sr;
        POINT *p;
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
	double x;

	infile = fopen(restart_name,"r");

        /* Initialize states at the interface */
        next_output_line_containing_string(infile,"Interface states:");
        next_point(intfc,NULL,NULL,NULL);
        while (next_point(intfc,&p,&hse,&hs))
        {
            FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl,(POINTER*)&sr);
            fscanf(infile,"%lf",&x);
            *sl = (STATE)x;
            fscanf(infile,"%lf",&x);
            *sr = (STATE)x;
        }

	FT_MakeGridIntfc(front);
        setDomain();

        /* Initialize states in the interior regions */

	next_output_line_containing_string(infile,"Interior states:");

	switch (dim)
	{
	case 1:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
		index = d_index1d(i,top_gmax);
	    	fscanf(infile,"%lf",&cell_center[index].state.T);
	    }
	    break;
	case 2:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    {
		index = d_index2d(i,j,top_gmax);
	    	fscanf(infile,"%lf",&cell_center[index].state.T);
	    }
	    break;
	case 3:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (k = 0; k <= top_gmax[2]; ++k)
	    {
		index = d_index3d(i,j,k,top_gmax);
	    	fscanf(infile,"%lf",&cell_center[index].state.T);
	    }
	    break;
	}
	fclose(infile);
	copyMeshStates();
}

void CARTESIAN::setBoundary()
{
	int i,j,k,index0,index1;
	INTERFACE *intfc = front->interf;

	switch (dim)
	{
	case 1:
	    if (rect_boundary_type(intfc,0,0) == NEUMANN_BOUNDARY)
	    {
		index0 = d_index1d(0,top_gmax);
		index1 = d_index1d(1,top_gmax);
	    	cell_center[index0].state.T  = cell_center[index1].state.T;
	    }
	    if (rect_boundary_type(intfc,0,1) == NEUMANN_BOUNDARY)
	    {
		index0 = d_index1d(top_gmax[0],top_gmax);
		index1 = d_index1d(top_gmax[0]-1,top_gmax);
	    	cell_center[index0].state.T  = cell_center[index1].state.T;
	    }
	    break;
	case 2:
	    if (rect_boundary_type(intfc,0,0) == NEUMANN_BOUNDARY)
	    {
		for (j = 0; j <= top_gmax[1]; ++j)
		{
		    index0 = d_index2d(0,j,top_gmax);
		    index1 = d_index2d(1,j,top_gmax);
	    	    cell_center[index0].state.T  = cell_center[index1].state.T;
		}
	    }
	    if (rect_boundary_type(intfc,0,1) == NEUMANN_BOUNDARY)
	    {
		for (j = 0; j <= top_gmax[1]; ++j)
		{
		    index0 = d_index2d(top_gmax[0],j,top_gmax);
		    index1 = d_index2d(top_gmax[0]-1,j,top_gmax);
	    	    cell_center[index0].state.T  = cell_center[index1].state.T;
		}
	    }
	    if (rect_boundary_type(intfc,1,0) == NEUMANN_BOUNDARY)
	    {
		for (i = 0; i <= top_gmax[0]; ++i)
		{
		    index0 = d_index2d(i,0,top_gmax);
		    index1 = d_index2d(i,1,top_gmax);
	    	    cell_center[index0].state.T  = cell_center[index1].state.T;
		}
	    }
	    if (rect_boundary_type(intfc,1,1) == NEUMANN_BOUNDARY)
	    {
		for (i = 0; i <= top_gmax[0]; ++i)
		{
		    index0 = d_index2d(i,top_gmax[1],top_gmax);
		    index1 = d_index2d(i,top_gmax[1]-1,top_gmax);
	    	    cell_center[index0].state.T  = cell_center[index1].state.T;
		}
	    }
	    break;
	case 3:
	    if (rect_boundary_type(intfc,0,0) == NEUMANN_BOUNDARY)
	    {
		for (j = 0; j <= top_gmax[1]; ++j)
		for (k = 0; k <= top_gmax[2]; ++k)
		{
		    index0 = d_index3d(0,j,k,top_gmax);
		    index1 = d_index3d(1,j,k,top_gmax);
	    	    cell_center[index0].state.T  = cell_center[index1].state.T;
		}
	    }
	    if (rect_boundary_type(intfc,0,1) == NEUMANN_BOUNDARY)
	    {
		for (j = 0; j <= top_gmax[1]; ++j)
		for (k = 0; k <= top_gmax[2]; ++k)
		{
		    index0 = d_index3d(top_gmax[0],j,k,top_gmax);
		    index1 = d_index3d(top_gmax[0]-1,j,k,top_gmax);
	    	    cell_center[index0].state.T  = cell_center[index1].state.T;
		}
	    }
	    if (rect_boundary_type(intfc,1,0) == NEUMANN_BOUNDARY)
	    {
		for (i = 0; i <= top_gmax[0]; ++i)
		for (k = 0; k <= top_gmax[2]; ++k)
		{
		    index0 = d_index3d(i,0,k,top_gmax);
		    index1 = d_index3d(i,1,k,top_gmax);
	    	    cell_center[index0].state.T  = cell_center[index1].state.T;
		}
	    }
	    if (rect_boundary_type(intfc,1,1) == NEUMANN_BOUNDARY)
	    {
		for (i = 0; i <= top_gmax[0]; ++i)
		for (k = 0; k <= top_gmax[2]; ++k)
		{
		    index0 = d_index3d(i,top_gmax[1],k,top_gmax);
		    index1 = d_index3d(i,top_gmax[1]-1,k,top_gmax);
	    	    cell_center[index0].state.T  = cell_center[index1].state.T;
		}
	    }
	    if (rect_boundary_type(intfc,2,0) == NEUMANN_BOUNDARY)
	    {
		for (i = 0; i <= top_gmax[0]; ++i)
		for (j = 0; j <= top_gmax[1]; ++j)
		{
		    index0 = d_index3d(i,j,0,top_gmax);
		    index1 = d_index3d(i,j,1,top_gmax);
	    	    cell_center[index0].state.T  = cell_center[index1].state.T;
		}
	    }
	    if (rect_boundary_type(intfc,2,1) == NEUMANN_BOUNDARY)
	    {
		for (i = 0; i <= top_gmax[0]; ++i)
		for (j = 0; j <= top_gmax[1]; ++j)
		{
		    index0 = d_index3d(i,j,top_gmax[2],top_gmax);
		    index1 = d_index3d(i,j,top_gmax[2]-1,top_gmax);
	    	    cell_center[index0].state.T  = cell_center[index1].state.T;
		}
	    }
	    break;
	}
}	/* end setBoundary */

void CARTESIAN::oneDimPlot(char *outname)
{

#if defined __GD__
	gdOneDimPlot(outname);
#endif /* defined __GD__ */
	xgraphOneDimPlot(outname);
}	/* end temperaturePlot */

void CARTESIAN::xgraphOneDimPlot(char *outname)
{
	int i,index;
	char filename[100];
	FILE *outfile;

	if (debugging("trace"))
	    printf("Entering xgraphTemp1()\n");
        sprintf(filename,"%s/tmp-xg.ts%s",outname,right_flush(front->step,7));
#if defined(__MPI__)
        sprintf(filename,"%s-nd%s",filename,right_flush(pp_mynode(),4));
#endif /* defined(__MPI__) */
        outfile = fopen(filename,"w");
	fprintf(outfile,"\"Solid temp at %6.3f\"\n",front->time);
	for (i = 0; i <= top_gmax[0]; ++i)
	{
	    index = d_index1d(i,top_gmax);
	    if (cell_center[index].comp == SOLID_COMP)
	    	fprintf(outfile,"%24.18g  %24.18g\n",
	    		cell_center[index].coords[0],
	    		cell_center[index].state.T);
	}
	fprintf(outfile,"\n\n");
	fprintf(outfile,"\"Liquid temp at %6.3f\"\n",front->time);
	for (i = 0; i <= top_gmax[0]; ++i)
	{
	    index = d_index1d(i,top_gmax);
	    if (cell_center[index].comp == LIQUID_COMP)
	    	fprintf(outfile,"%24.18g  %24.18g\n",
	    		cell_center[index].coords[0],
	    		cell_center[index].state.T);
	}
	fclose(outfile);
	if (debugging("trace"))
	    printf("Leaving xgraphTemp1()\n");
}	/* end xgraphOneDimPlot */


void CARTESIAN::computeAdvectionExplicit(COMPONENT sub_comp)
{
	int i,j,k,l,m,ic,icn,icoords[MAXD],nc;
	int gmin[MAXD],ipn[MAXD];
	int index0;
	double coords[MAXD],crx_coords[MAXD];
	double temperature,temperature_nb,dgrad[MAXD];
	double coef;
	COMPONENT comp;
	boolean fr_crx_grid_seg;
	const GRID_DIRECTION dir[3][2] = 
		{{WEST,EAST},{SOUTH,NORTH},{LOWER,UPPER}};

	start_clock("computeAdvectionExplicit");

	coef = eqn_params->D*m_dt;

	for (i = 0; i < dim; ++i) gmin[i] = 0;

	switch (dim)
	{
        case 1:
            for (i = imin; i <= imax; ++i)
            {
                icoords[0] = i;
                ic = d_index1d(i,top_gmax);
                comp = top_comp[ic];
                if (comp != sub_comp)
                     continue;
                array[ic] = temperature = cell_center[ic].state.T;
                for (l = 0; l < dim; ++l)
                {
                    dgrad[l] = 0.0;
                    for (m = 0; m < 2; ++m)
                    {
                        fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                icoords,dir[l][m],comp,temperature_of_state,
                                &temperature_nb,crx_coords);
                        if (!fr_crx_grid_seg)
                        {
                             next_ip_in_dir(icoords,dir[l][m],ipn,gmin,
							top_gmax);
                             icn = d_index1d(ipn[0],top_gmax);
			     temperature_nb = cell_center[icn].state.T;
                        }
                        dgrad[l] += (temperature_nb - temperature)/top_h[l];
                    }
                    array[ic] += coef*dgrad[l]/top_h[l];
                }
            }
            break;
	case 2:
	    for (i = imin; i <= imax; ++i)
	    for (j = jmin; j <= jmax; ++j)
	    {
	    	icoords[0] = i;
	    	icoords[1] = j;
	    	ic = d_index2d(i,j,top_gmax);
	    	comp = top_comp[ic];
	    	if (comp != sub_comp) 
	    	    continue;
                array[ic] = temperature = cell_center[ic].state.T;
		for (l = 0; l < dim; ++l)
		{
	            dgrad[l] = 0.0;
                    for (m = 0; m < 2; ++m)
                    {
                        fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                icoords,dir[l][m],comp,temperature_of_state,
                                &temperature_nb,crx_coords);
                        if (!fr_crx_grid_seg)
                        {
                            next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                            icn = d_index2d(ipn[0],ipn[1],top_gmax);
			    temperature_nb = cell_center[icn].state.T;
                        }
                        dgrad[l] += (temperature_nb - temperature)/top_h[l];
                    }
		    array[ic] += coef*dgrad[l]/top_h[l];
		}
	    }
	    break;
	case 3:
	    for (i = imin; i <= imax; ++i)
	    for (j = jmin; j <= jmax; ++j)
	    for (k = kmin; k <= kmax; ++k)
	    {
	    	icoords[0] = i;
	    	icoords[1] = j;
	    	icoords[2] = k;
	    	ic = d_index3d(i,j,k,top_gmax);
	    	comp = top_comp[ic];
	    	if (comp != sub_comp) 
	    	    continue;
                array[ic] = temperature = cell_center[ic].state.T;
		for (l = 0; l < dim; ++l)
		{
	            dgrad[l] = 0.0;
                    for (m = 0; m < 2; ++m)
                    {
                        fr_crx_grid_seg = FT_StateVarAtGridCrossing(front,
                                icoords,dir[l][m],comp,temperature_of_state,
                                &temperature_nb,crx_coords);
                        if (!fr_crx_grid_seg)
                        {
                            next_ip_in_dir(icoords,dir[l][m],ipn,gmin,top_gmax);
                            icn = d_index3d(ipn[0],ipn[1],ipn[2],top_gmax);
                	    temperature_nb = cell_center[icn].state.T;
                        }
                        dgrad[l] += (temperature_nb - temperature)/top_h[l];
                    }
		    array[ic] += coef*dgrad[l]/top_h[l];
		}
	    }
	    break;
	}
	scatMeshArray();
	switch (dim)
	{
        case 1:
            for (i = 0; i <= top_gmax[0]; ++i)
            {
		ic = d_index1d(i,top_gmax);
	    	comp = top_comp[ic];
		if (comp == sub_comp)
		    cell_center[ic].state.T = array[ic];
	    }
	    break;
        case 2:
            for (i = 0; i <= top_gmax[0]; ++i)
            for (j = 0; j <= top_gmax[1]; ++j)
            {
		ic = d_index2d(i,j,top_gmax);
	    	comp = top_comp[ic];
		if (comp == sub_comp)
		    cell_center[ic].state.T = array[ic];
	    }
	    break;
        case 3:
            for (i = 0; i <= top_gmax[0]; ++i)
            for (j = 0; j <= top_gmax[1]; ++j)
            for (k = 0; k <= top_gmax[2]; ++k)
            {
		ic = d_index3d(i,j,k,top_gmax);
	    	comp = top_comp[ic];
		if (comp == sub_comp)
		    cell_center[ic].state.T = array[ic];
	    }
	    break;
	}

	stop_clock("computeAdvectionExplicit");
}	/* computeAdvectionExplicit */

#if defined __GD__
void CARTESIAN::gdOneDimPlot(
	char *out_name)
{
	int step = front->step;
	INTERFACE *intfc = front->interf;
	int i,index0;
        double *x,*c;
        char movie_caption[100];
        char time_label[100];
        char gd_name[200];
        double xmin,xmax,cmin,cmax,height;
	static boolean first = YES;
	POINT *p,*pts[100];
        HYPER_SURF *hs;
        HYPER_SURF_ELEMENT *hse;
        STATE *sl[100],*sr[100];
	int np,count;
	boolean tracked_interior_point;

	if (debugging("trace"))
	    printf("Entering gdOneDimPlot()\n");

	FT_VectorMemoryAlloc((POINTER*)&x,top_gmax[0]+1,sizeof(double));
	FT_VectorMemoryAlloc((POINTER*)&c,top_gmax[0]+1,sizeof(double));

	if (first)
	{
	    first = NO;
	    cmin = HUGE;	cmax = -HUGE;
	    for (i = 1; i < top_gmax[0]; ++i)
	    {
	    	index0 = d_index1d(i,top_gmax);
	    	c[i] = cell_center[index0].state.T;
		if (cmin > c[i]) cmin = c[i];
            	if (cmax < c[i]) cmax = c[i];
	    }
	    height = cmax - cmin;
	    xmin = top_L[0];	xmax = top_U[0];
	    cmin -= 0.15*height;    cmax += 0.30*height;
	    sprintf(movie_caption,"T vs. x");
            sprintf(gd_name,"%s/temp-gd.gif",out_name);
            gd_initplot(gd_name,movie_caption,xmin,xmax,cmin,cmax,3);
	}

	np = 0;
	next_point(intfc,NULL,NULL,NULL);
	while (next_point(intfc,&p,&hse,&hs))
	{
	    pts[np] = p;
	    FT_GetStatesAtPoint(p,hse,hs,(POINTER*)&sl[np],(POINTER*)&sr[np]);
	    np++;
	}

	count = np = 0;
	for (i = 1; i < top_gmax[0]; ++i)
	{
	    x[count] = top_L[0] + i*top_h[0];
	    if (i == 1 && x[count] > Coords(pts[np])[0]) np++;
	    index0 = d_index1d(i,top_gmax);
	    c[count] = cell_center[index0].state.T;
	    if (x[count] > Coords(pts[np])[0]) 
	    {
		x[count] = Coords(pts[np])[0];
	    	c[count++] = (double)*sl[np];
		gd_plotdata(count,x,c);
		count = 0;
		x[count] = Coords(pts[np])[0];
		c[count++] = (double)*sr[np];
		np++;
	    }
	    else if (i == top_gmax[0] - 1)
	    {
		count++;
		gd_plotdata(count,x,c);
	    }
	    else
	    	count++;
	}

	sprintf(time_label,"Time = %6.3f",front->time);
	gd_plotframe(time_label);
	if (debugging("trace"))
	    printf("Leaving gdOneDimPlot()\n");

}	/* end gdOneDimPlot */
#endif /* defined __GD__ */


void CARTESIAN::setDomain()
{
	static boolean first = YES;
	INTERFACE *grid_intfc;
	Table *T;
	int i;

	if (debugging("trace")) printf("Passed FT_MakeGridIntfc()\n");

	grid_intfc = front->grid_intfc;
	top_grid = &topological_grid(grid_intfc);
	lbuf = front->rect_grid->lbuf;
	ubuf = front->rect_grid->ubuf;
	top_gmax = top_grid->gmax;
	top_L = top_grid->L;
	top_U = top_grid->U;
	top_h = top_grid->h;
	dim = grid_intfc->dim;
	T = table_of_interface(grid_intfc);
	top_comp = T->components;
	eqn_params = (PARAMS*)front->extra2;
	hmin = top_h[0];
	for (i = 1; i < dim; ++i)
	    if (hmin > top_h[i]) hmin = top_h[i];

	switch (dim)
	{
	case 1:
            if (first)
            {
		FT_ScalarMemoryAlloc((POINTER*)&field,sizeof(PHASE_FIELD));
                FT_VectorMemoryAlloc((POINTER*)&array,top_gmax[0]+1,FLOAT);
                FT_VectorMemoryAlloc((POINTER*)&field->temperature,top_gmax[0]+1,FLOAT);
                first = NO;
            }	
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    eqn_params->field = field;
	    break;
	case 2:
	    if (first)
	    {
		FT_ScalarMemoryAlloc((POINTER*)&field,sizeof(PHASE_FIELD));
	    	FT_VectorMemoryAlloc((POINTER*)&array,(top_gmax[0]+1)*(top_gmax[1]+1),FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->temperature,
			(top_gmax[0]+1)*(top_gmax[1]+1),FLOAT);
	    	first = NO;
	    }
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
	    eqn_params->field = field;
	    break;
	case 3:
	    if (first)
	    {
		FT_ScalarMemoryAlloc((POINTER*)&field,sizeof(PHASE_FIELD));
	    	FT_VectorMemoryAlloc((POINTER*)&array,
			(top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1),FLOAT);
	    	FT_VectorMemoryAlloc((POINTER*)&field->temperature,
			(top_gmax[0]+1)*(top_gmax[1]+1)*(top_gmax[2]+1),FLOAT);
	    	first = NO;
	    }
	    imin = (lbuf[0] == 0) ? 1 : lbuf[0];
	    jmin = (lbuf[1] == 0) ? 1 : lbuf[1];
	    kmin = (lbuf[2] == 0) ? 1 : lbuf[2];
	    imax = (ubuf[0] == 0) ? top_gmax[0] - 1 : top_gmax[0] - ubuf[0];
	    jmax = (ubuf[1] == 0) ? top_gmax[1] - 1 : top_gmax[1] - ubuf[1];
	    kmax = (ubuf[2] == 0) ? top_gmax[2] - 1 : top_gmax[2] - ubuf[2];
	    eqn_params->field = field;
	    break;
	}
}	/* end setDomain */

void CARTESIAN::copyMeshStates()
{
	int i,j,k,index;
	double *temperature = eqn_params->field->temperature;

	switch (dim)
	{
	case 1:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    {
	    	index  = d_index1d(i,top_gmax);
	    	temperature[index] = cell_center[index].state.T;
	    }
	    break;
	case 2:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    {
	    	index  = d_index2d(i,j,top_gmax);
	    	temperature[index] = cell_center[index].state.T;
	    }
	    break;
	case 3:
	    for (i = 0; i <= top_gmax[0]; ++i)
	    for (j = 0; j <= top_gmax[1]; ++j)
	    for (k = 0; k <= top_gmax[2]; ++k)
	    {
	    	index  = d_index3d(i,j,k,top_gmax);
	    	temperature[index] = cell_center[index].state.T;
	    }
	    break;
	}
}	/* end copySolute */
