#include "callbackMtoS.h"

callbackMtoS::callbackMtoS(GRBModel *m, Params *params){
	callback_params = params;
	model = m;
	cout<<"In callbackMtoS"<<endl;
}

void callbackMtoS::callback(){
	// cout<<"Im callback"<<endl;
	if(sol_repeat_counter<100){
		Matrix<double,3,1> T_init_ICP;
		Matrix<double,3,3> temp_R;
		Matrix<double,3,3> R_final;
		Matrix<double,1,3> T_final;

		if(where==GRB_CB_MIPNODE){
			int numvars = model->get(GRB_IntAttr_NumVars);
		    GRBVar *vars = model->getVars();
			double *sol_at_node = getNodeRel(vars, numvars);
			int l = 0;
			for(int i=4*(*(callback_params->Ns_sampled)); i<4*(*(callback_params->Ns_sampled))+3; i++){
				T_init_ICP(l,0) = sol_at_node[i];
				temp_R(0,l) = sol_at_node[i+3];
				temp_R(1,l) = sol_at_node[i+6];
				temp_R(2,l) = sol_at_node[i+9];
				l = l+1;
			}
			JacobiSVD<MatrixXd> svd(temp_R, ComputeThinU | ComputeThinV);	// Computes U,Vt.T,diagonal of Sigma Matrices.
			MatrixXd R_init_ICP = svd.matrixU()*svd.matrixV().transpose();
			if(R_init_ICP.determinant()<0){
				MatrixXd temp = svd.matrixV();
				temp(0,2)=-1*temp(0,2); temp(1,2)=-1*temp(1,2); temp(2,2)=-1*temp(2,2);
				R_init_ICP = svd.matrixU()*temp.transpose();
			}

			Matrix<double,4,4> temp_transf_mat = MatrixXd::Zero(4,4);
			temp_transf_mat.block(0,0,3,3)=R_init_ICP;
			temp_transf_mat.block(0,3,3,1)=T_init_ICP;
			temp_transf_mat(3,3) = 1;
			
			// Take inverse of the temp_transf_mat to find rotation and translation.
			R_init_ICP = temp_transf_mat.inverse().block(0,0,3,3);
			T_init_ICP = temp_transf_mat.inverse().block(0,3,3,1);

			// transformed_for_ICP = R_init_ICP * S[:,num_sampled_sens_points_ICP] + T_init_ICP
			MatrixXd transformed_for_ICP = (R_init_ICP*((*(callback_params->S)).block(0,0,3,*(callback_params->num_sampled_sens_points_ICP)))).colwise()+T_init_ICP;

			// Run ICP on the transformed sensor data and model points.
			MatrixXd transf_ICP_output = icp_test(callback_params->V, &transformed_for_ICP, callback_params->M, callback_params->tree_M, callback_params->M_sampled, callback_params->tree_M_sampled, callback_params->SigmaS, callback_params->F, callback_params->points_per_face, callback_params->ICP_triangle_proj_switch_callback, callback_params->ICP_or_GICP_switch_callback);

			// Transformation Matrix from Sensor to Model: (R_ICP_output,T_ICP_output)*(R_init_ICP,T_init_ICP)
			MatrixXd total_transf_ICP_S_to_M = transf_ICP_output*temp_transf_mat;

			// Input to global optimizer.
			MatrixXd total_transf_ICP_M_to_S = total_transf_ICP_S_to_M.inverse();

			Matrix<double,3,1> temp_T_S_to_M = total_transf_ICP_S_to_M.block(0,3,3,1);	// Translation part from Sensor to Model.
			// transformed_S_sampled_after_ICP = total_transf_ICP_S_to_M[0:3,0:3]*S[:,0:Ns_sampled]+temp_T_S_to_M
			MatrixXd transformed_S_sampled_after_ICP = ((total_transf_ICP_S_to_M.block(0,0,3,3)) * ((*(callback_params->S)).block(0,0,3,*(callback_params->Ns_sampled)))).colwise()+temp_T_S_to_M;
			
			OptVariables opt_vars_callback;			// find optimization variables in callback.
			// find all documentation in the helper function.
			find_Cb(&opt_vars_callback.Cb, &transformed_S_sampled_after_ICP, callback_params->tree_M_sampled, (*(callback_params->M_global)).cols());
			find_lam(&opt_vars_callback.lam, &total_transf_ICP_S_to_M, callback_params->num_partitions_SOS2);	// total_transf_ICP_S_to_M because this function takes inverse of given matrix and then find the lambda.
			find_w(&opt_vars_callback.w, &opt_vars_callback.lam, callback_params->num_partitions_SOS2);

			MatrixXd S_alpha_ip = (*(callback_params->S)).block(0,0,3,*(callback_params->Ns_sampled));		// Input Sensor data to find_alpha function.
			find_alpha(&opt_vars_callback.alpha, &total_transf_ICP_S_to_M, &S_alpha_ip, callback_params->M_global, &opt_vars_callback.Cb, callback_params->B);		// total_transf_ICP_S_to_M because this function takes inverse of given matrix and then find the lambda.
			find_phi(&opt_vars_callback.phi, &opt_vars_callback.alpha);

			cout<<"Phi After ICP in callback: "<<(opt_vars_callback.phi.sum())/(*(callback_params->Ns_sampled))<<endl;

			MatrixXd T_after_ICP_M_to_S = total_transf_ICP_M_to_S.block(0,3,3,1);	// find translation part.
			MatrixXd R_after_ICP_M_to_S = total_transf_ICP_M_to_S.block(0,0,3,3);	// find rotation part.

			double list_T_after_ICP_M_to_S[T_after_ICP_M_to_S.cols()*T_after_ICP_M_to_S.rows()];			// list to store translations
			double list_R_after_ICP_M_to_S[R_after_ICP_M_to_S.cols()*R_after_ICP_M_to_S.rows()];			// list to store rotation
			double list_Cb_sampled_after_ICP[opt_vars_callback.Cb.cols()*opt_vars_callback.Cb.rows()];		// list to store correspondences
			double list_w_after_ICP[opt_vars_callback.w.cols()*opt_vars_callback.w.rows()];					// list to store w
			double list_phi_after_ICP[opt_vars_callback.phi.cols()*opt_vars_callback.phi.rows()];			// list to store phi
			double list_alpha_after_ICP[opt_vars_callback.alpha.cols()*opt_vars_callback.alpha.rows()];		// list to store alpha
			double list_lam_after_ICP[opt_vars_callback.lam[0].cols()*opt_vars_callback.lam[0].rows()*opt_vars_callback.lam.size()];	// list to store lam

			// Flatten the matrices to store in the list. (Useful for setSolution)
			flatten_matrix(&T_after_ICP_M_to_S, list_T_after_ICP_M_to_S);
			flatten_matrix(&R_after_ICP_M_to_S, list_R_after_ICP_M_to_S);
			flatten_matrix(&opt_vars_callback.Cb, list_Cb_sampled_after_ICP);
			flatten_w(&opt_vars_callback.w, list_w_after_ICP);
			flatten_matrix(&opt_vars_callback.phi, list_phi_after_ICP);
			flatten_matrix(&opt_vars_callback.alpha, list_alpha_after_ICP);
			flatten_vector(&opt_vars_callback.lam, list_lam_after_ICP);

			// setSolution for T.
			int start_count = 4*(*(callback_params->Ns_sampled));
			GRBVar T_vars[T_after_ICP_M_to_S.cols()*T_after_ICP_M_to_S.rows()];		// Define gurobi variables array.
			for(int i=0; i<T_after_ICP_M_to_S.cols()*T_after_ICP_M_to_S.rows(); i++){
				T_vars[i]=vars[start_count+i];											// Store each variable in variable array.
			}
			setSolution(T_vars, list_T_after_ICP_M_to_S, (T_after_ICP_M_to_S.cols()*T_after_ICP_M_to_S.rows()));	// set the values of variables found in callback.

			// setSolution for R.
			start_count = 4*(*(callback_params->Ns_sampled))+3;
			GRBVar R_vars[R_after_ICP_M_to_S.cols()*R_after_ICP_M_to_S.rows()];
			for(int i=0; i<R_after_ICP_M_to_S.cols()*R_after_ICP_M_to_S.rows(); i++){
				R_vars[i]=vars[start_count+i];
			}
			setSolution(R_vars, list_R_after_ICP_M_to_S, (R_after_ICP_M_to_S.cols()*R_after_ICP_M_to_S.rows()));

			// setSolution for Cb.
			start_count = 4*(*(callback_params->Ns_sampled))+3+9+9;
			GRBVar Cb_vars[opt_vars_callback.Cb.cols()*opt_vars_callback.Cb.rows()];
			for(int i=0; i<opt_vars_callback.Cb.cols()*opt_vars_callback.Cb.rows(); i++){
				Cb_vars[i]=vars[start_count+i];
			}
			setSolution(Cb_vars, list_Cb_sampled_after_ICP, opt_vars_callback.Cb.cols()*opt_vars_callback.Cb.rows());

			// setSolution for alpha.
			start_count = 0;
			GRBVar alpha_vars[opt_vars_callback.alpha.cols()*opt_vars_callback.alpha.rows()];
			for(int i=0; i<opt_vars_callback.alpha.cols()*opt_vars_callback.alpha.rows(); i++){
				alpha_vars[i]=vars[start_count+i];
			}
			setSolution(alpha_vars, list_alpha_after_ICP, opt_vars_callback.alpha.cols()*opt_vars_callback.alpha.rows());

			// setSolution for phi.
			start_count = 3*(*(callback_params->Ns_sampled));
			GRBVar phi_vars[opt_vars_callback.phi.cols()*opt_vars_callback.phi.rows()];
			for(int i=0; i<opt_vars_callback.phi.cols()*opt_vars_callback.phi.rows(); i++){
				phi_vars[i]=vars[start_count+i];
			}
			setSolution(phi_vars, list_phi_after_ICP, opt_vars_callback.phi.cols()*opt_vars_callback.phi.rows());

			// setSolution for lam.
			start_count = 4*(*(callback_params->Ns_sampled)) +3 +9 +9 +(*(callback_params->Ns_sampled)) * ((*(callback_params->M_global)).cols());
			GRBVar lam_vars[opt_vars_callback.lam[0].cols()*opt_vars_callback.lam[0].rows()*opt_vars_callback.lam.size()];
			for(int i=0; i<opt_vars_callback.lam[0].cols()*opt_vars_callback.lam[0].rows()*(int)(opt_vars_callback.lam.size()); i++){
				lam_vars[i]=vars[start_count+i];
			}
			setSolution(lam_vars, list_lam_after_ICP, opt_vars_callback.lam[0].cols()*opt_vars_callback.lam[0].rows()*opt_vars_callback.lam.size());

			double objval = useSolution();			// useSolution method.

			if(objval == 1e+100){
				sol_repeat_counter = sol_repeat_counter+1;
			}
			else{
				sol_repeat_counter = 0;
			}

			cout<<"Objective Value: "<<objval<<endl;
			cout<<"Sol Repeat Counter: "<<sol_repeat_counter<<endl;
		}
	}
}