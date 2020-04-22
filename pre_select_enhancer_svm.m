classdef pre_select_enhancer_svm < handle
	
	properties
		pos_sample_inputs;
		pos_sample_outputs;
		neg_sample_inputs;
		neg_sample_outputs;
        
        testing_outputs;
        testing_predictions;
		
		kmer_length;
		
		positive_file;
		negative_file;
        
        features_name;
        weight;
        
                
        pos_use;
        neg_use;
	end
	
	methods
		
		function a = pre_select_enhancer_svm(positive_file,negative_file,kmer_length,pos_use,neg_use)
			a.kmer_length = kmer_length;
			
			a.positive_file = positive_file;
			a.negative_file = negative_file;
			            
            a.pos_use = pos_use;
            a.neg_use = neg_use;
            
			%Load data
			a.get_data();
			
            a.single_fold();
			
			%a.plot_roc()
        end
		
        function single_fold(a)
            pos_test_set_outputs = a.pos_sample_outputs;
            neg_test_set_outputs = a.neg_sample_outputs;
            pos_train_set_outputs = a.pos_sample_outputs;
            neg_train_set_outputs = a.neg_sample_outputs;

            testing_inputs = [a.pos_sample_inputs;a.neg_sample_inputs];
            a.testing_outputs = [pos_test_set_outputs;neg_test_set_outputs];

            training_inputs = [a.pos_sample_inputs;a.neg_sample_inputs];
            training_outputs = [pos_train_set_outputs;neg_train_set_outputs];

            disp(['Size of training feartures vectors: ',num2str(size(training_inputs,1))]);

            %Train the svm
            svm = a.train(training_inputs,training_outputs);

            %Test the svm
            [a.testing_predictions,a.weight] = a.test(svm,testing_inputs);
        end
		
		%Plot an error against ensemble size graph in a separate window.
		function [ROC_AUC,PR_AUC,feature_weight] = plot_roc(a,auc_only)
			n = length(a.testing_predictions);
            if(~auc_only)
                hold all;
                line_style = {'-','--','-.',':','-',':'};
                %color = {'k','k','k','k','k','r'};
            end
            
			a_weight = a.weight;
                
            % fit probabilities for scores
            [ROC_X,ROC_Y,~,R_AUC] = perfcurve(a.testing_outputs,a.testing_predictions,1);
            ROC_AUC = R_AUC;

            [PR_X,PR_Y,~,P_AUC] = perfcurve(a.testing_outputs,a.testing_predictions,1,'xCrit','reca','yCrit','prec');
            PR_AUC = P_AUC;

            if(~auc_only)
                % plot only if needed
                p = plot(ROC_X,ROC_Y);

                % Set DisplayNames for the lines for use by the legend
                %disp(['Fold ',num2str(i),'/',num2str(n),' AUC:',num2str(AUC)]);
                set(p, {'LineStyle'}, {line_style{1}});
                %set(p, {'Color'},{color{i}});
                set(p(1),'Displayname',['Fold ',num2str(1),'/',num2str(n),' ROC AUC:',num2str(R_AUC)]);

                % plot only if needed
                p = plot(PR_X,PR_Y);

                % Set DisplayNames for the lines for use by the legend
                %disp(['Fold ',num2str(i),'/',num2str(n),' AUC:',num2str(AUC)]);
                set(p, {'LineStyle'}, {line_style{1}});
                %set(p, {'Color'},{color{i}});
                set(p(1),'Displayname',['Fold ',num2str(1),'/',num2str(n),' PR AUC:',num2str(P_AUC)]);
            end
			
            feature_weight = [a.features_name',num2cell(a_weight')];

            feature_weight = sortrows(feature_weight,-2);
			
			disp(['The average ROC AUC is ',num2str(ROC_AUC)]);
            
            disp(['The average PR AUC is ',num2str(PR_AUC)]);
			
            if(~auc_only)
                % Center a legend at the top of the graph
                legend('Location','SouthEast');
                xlabel('False positive rate and Recall'); ylabel('True positive rate and Precision');
                title('ROC and PR');

                hold off;
            end
        end
		
		%Train SVM
        function svm =  train(a,training_inputs,training_outputs)
			%a.svm_struct{i} = svmtrain(training_inputs,training_outputs,'boxconstraint',0.00001,'kernel_function','linear','options',statset('MaxIter',500000));
			%a.svm_struct{i} = svmtrain(training_inputs,training_outputs,'autoscale',false,'boxconstraint',0.0001,'kernel_function','linear','method','LS');
            
            %a.svm_struct{i} = svmtrain(training_inputs,training_outputs,'autoscale',false,'boxconstraint',100,'kernel_function','linear','method','LS');
            svm = svmtrain(training_inputs,training_outputs,'autoscale',false,'boxconstraint',100,'kernel_function','linear','method','SMO','kernelcachelimit',50000,'options',statset('MaxIter',100000));
            %svm = svmtrain(training_inputs,training_outputs,'autoscale',false,'boxconstraint',0.1,'kernel_function','linear','method','SMO','kernelcachelimit',50000,'options',statset('MaxIter',100000));
		end
		
		%Test SVM
		function [test_outputs,weight]=test(a,svm,test_inputs)
			%test_outputs = svmclassify(a.svm_struct{i},test_inputs);
			
            %{
            shift = svm.ScaleData.shift;
			scale = svm.ScaleData.scaleFactor;
			test_inputs = bsxfun(@plus,test_inputs,shift);
			test_inputs = bsxfun(@times,test_inputs,scale);
            %}
            
			sv = svm.SupportVectors;
			alphaHat = svm.Alpha;
			bias = svm.Bias;
			kfun = svm.KernelFunction;
			kfunargs = svm.KernelFunctionArgs;
			test_outputs = kfun(sv,test_inputs,kfunargs{:})'*alphaHat(:) + bias;
            %test_outputs = kfun(sv,test_inputs,kfunargs{:})'*alphaHat(:);
			test_outputs = -test_outputs; % flip the sign to get the score for the +1 class
            
            weight  = -(alphaHat' * sv);
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
