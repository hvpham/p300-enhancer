classdef pre_select_enhancer < handle
	
	properties
		pos_sample_inputs;
		pos_sample_outputs;
		neg_sample_inputs;
		neg_sample_outputs;

        kmer_length;
        method;
        no_preselected;
		
		positive_file;
		negative_file;

        pos_use;
        neg_use;
        
        features_name;
        
        selected_features;
	end
	
	methods
		
		function a = pre_select_enhancer(positive_file,negative_file,kmer_length,method,no_preselected,pos_use,neg_use)
			a.kmer_length = kmer_length;
            a.method = method;
            a.no_preselected = no_preselected;
			
			a.positive_file = positive_file;
			a.negative_file = negative_file;
			            
            a.pos_use = pos_use;
            a.neg_use = neg_use;
            
			%Load data
			a.get_data();
            
            %preselect
            a.pre_selected_feature()
            
        end
        
        function pre_selected_feature(a)
            inputs = [a.pos_sample_inputs;a.neg_sample_inputs];
            outputs = (-[a.pos_sample_outputs;a.neg_sample_outputs]/2)+1.5;
            
            if(strcmp(a.method, 'mRMR'))
                out = a.pre_selected_feature_mRMR(inputs,outputs);
            elseif(strcmp(a.method, 'tTest'))
                out = fsTtest(inputs,outputs);
            elseif(strcmp(a.method, 'infoGain'))
                out = fsInfoGain(inputs,outputs);
            elseif(strcmp(a.method, 'SBMLR'))
                out = fsSBMLR(inputs,outputs);
            elseif(strcmp(a.method, 'Fisher'))
                out = fsFisher(inputs,outputs);
            else
                error(['Do not support preselect: ',a.method]);
            end
            
            a.selected_features = a.features_name(out.fList);
        end
        
        function [out]=pre_selected_feature_mRMR(a,inputs,outputs)
            parm.k = a.no_preselected;
            parm.pool = 1000;
            parm.type = 1;
            out = fsMRMR(inputs,outputs,parm);
        end
        
        %Generate some data from the required function.
		function get_data(a)
			[a.features_name,a.pos_sample_inputs,a.pos_sample_outputs,a.neg_sample_inputs,a.neg_sample_outputs] = data_preproc(a.kmer_length,a.positive_file,a.negative_file);
			
            %use only the portion allowed
            a.pos_sample_inputs = a.pos_sample_inputs(a.pos_use,:);
            a.pos_sample_outputs = a.pos_sample_outputs(a.pos_use,:);
            a.neg_sample_inputs = a.neg_sample_inputs(a.neg_use,:);
            a.neg_sample_outputs = a.neg_sample_outputs(a.neg_use,:);
            
			%normalize
			for i=1:size(a.pos_sample_inputs,1)
				a.pos_sample_inputs(i,:) = a.pos_sample_inputs(i,:)/sum(a.pos_sample_inputs(i,:));
			end
			
			for i=1:size(a.neg_sample_inputs,1)
				a.neg_sample_inputs(i,:) = a.neg_sample_inputs(i,:)/sum(a.neg_sample_inputs(i,:));
            end
		end
		
	end
	
end
