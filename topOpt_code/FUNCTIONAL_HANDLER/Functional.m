
%--------------------------------------------------------------------------
% CLASS CREATED TO HANDLE FUNCTIONALS
%--------------------------------------------------------------------------

classdef Functional
    %----------------------------------------------------------------------
    %% PROPERTIES
    %----------------------------------------------------------------------
    properties
        % general properties
        problem_name = "";
        test = "";
        name = "";
        custom_name;
        mu;

        %functional value properties
        n_it;
        functional;
        f0Init;
        f0Init_out_box;
        f0Init_parts;
        func;
        func_out_box;
        func_alpha;
        func_grad;
        func_vort;
        func_p;
        valid;
        changes;
        no_w_functional;
        no_w_f_alpha;
        no_w_f_grad;
        no_w_f_vort;
        no_w_f_p;

        % matrix functional values
        func_values;
    end
    % END PROPERTIES
    %----------------------------------------------------------------------

    %----------------------------------------------------------------------
    %% PUBLIC METHODS 
    %----------------------------------------------------------------------
    methods (Access = public)
        %-------------------------------------------
        % CONSTRUCTOR
        function self = Functional(pb_name, pb_test, pb_mu, pb_custom) 
            self.problem_name = pb_name;
            self.test = pb_test;
            self.mu = pb_mu;
            % start_path = "test_res/";
            % end_path = "";
            start_path = "../results/";
            end_path = "/Matlab_interface";
            
            self.name = start_path + pb_name + "/" + pb_test + end_path;
            self.custom_name = pb_custom;
            self = eval(self); % evaluate functional
        end
        %-------------------------------------------

        %-------------------------------------------
        % PRINT
        function new_nfig = print(self, nfig, what_to_print)
            % print starting from nfig+1
            % what_to_print = [print_func, print_changes, print_no_w_func]
            new_nfig = nfig;
            if (what_to_print(1) == 1)
                new_nfig = new_nfig+1;
                print_func(self, new_nfig);
            end
            if (what_to_print(2) == 1)
                new_nfig = new_nfig+1;
                print_changes(self, new_nfig);
            end
            if (what_to_print(3) == 1)
                new_nfig = new_nfig+1;
                print_no_weight_func(self, new_nfig);
            end
        end
        %-------------------------------------------

        %-------------------------------------------
        % PLOT FUNCTIONAL VALUES
        function print_func(self, nfig)
            it = 1:self.n_it;
            figure(nfig);
            hold on;
            plot(it, self.func, "-", "Color", "#000000", "LineWidth", 2);
            plot(it, self.func_out_box,"-", "Color", "#A2142F", "LineWidth", 2);
            plot(it, self.func_alpha,"-", "Color", "#D95319", "LineWidth", 2);
            plot(it, self.func_grad,"-", "Color", "#0072BD", "LineWidth", 2);
            plot(it, self.func_vort,"-", "Color", "#77AC30", "LineWidth", 2);
            plot(it, self.func_p, "-", "Color", "#EDB120", "LineWidth", 2);
            for i = 1 : floor(length(self.func) / 20) : (length(self.func))
                if self.valid(i) == 1
                    plot(it(i), self.func(i), 'bo', "LineWidth", 1.5);
                else
                    plot(it(i), self.func(i), 'ro', "LineWidth", 1.5);
                end
            end
            legend("Jtot", "Jout", "Jalpha", "Jgrad", "Jvort", "Jp");
            title("FUNCTIONAL");
            % xscale("log");
            hold off;
        end
        %-------------------------------------------
        function print_specific_func(self, nfig, what_to_print, abs_value, print_changes, color, max_it)
            it = 1:self.n_it;
            figure(nfig);
            hold on;
            func_case = Compare_Functs.get_func_case(what_to_print);
            val_to_print = self.func_values(:,func_case);
            [max_val, max_id] = max(val_to_print(4:end));
            if (abs_value == 1)
                val_to_print = val_to_print * self.f0Init;
            end
            
            if (max_it == 0)
                real_max_it = self.n_it;
            else
                real_max_it = max_it;
            end

            print_it = it;
            print_val = val_to_print;
            for i_it=self.n_it+1:real_max_it
                print_it(i_it) = i_it;
                print_val(i_it) = val_to_print(self.n_it);
            end

            plot(print_it(1:real_max_it), print_val(1:real_max_it), "-", "Color", color, "LineWidth", 6, 'DisplayName', self.custom_name);

            if (print_changes == 1)
                % for i = 1 : floor(real_max_it/20) : (real_max_it)
                for i = 1 : (max_id+4)
                % for i = 1 : (14)
                    if self.valid(i) == 1
                        plot(print_it(i), print_val(i), 'ko', "LineWidth", 1.5, "MarkerSize", 10, 'MarkerFaceColor',[0.8 0.8 0.8]);
                    else
                        plot(print_it(i), print_val(i), 'k^', "LineWidth", 1.5, "MarkerSize", 10, 'MarkerFaceColor',[0.8 0.8 0.8]);
                    end
                end
            end

            hold off;
        end
        %-------------------------------------------

        %-------------------------------------------
        % PLOT NO WEIGHTS FUNCTIONAL VALUES
        function print_no_weight_func(self, nfig)
            it = 1:self.n_it;
            figure(nfig);
            hold on;
            plot(it, self.no_w_f_alpha,"-", "Color", "#D95319", "LineWidth", 2);
            plot(it, self.no_w_f_grad,"-", "Color", "#0072BD", "LineWidth", 2);
            plot(it, self.no_w_f_vort,"-", "Color", "#77AC30", "LineWidth", 2);
            plot(it, self.no_w_f_p, "-", "Color", "#EDB120", "LineWidth", 2);
            legend("Jalpha", "Jgrad", "Jvort", "Jp");
            title("NO W FUNCTIONAL");
            hold off;
        end
        %-------------------------------------------

        %-------------------------------------------
        % PLOT FUNCTIONAL CHANGES
        function print_changes(self, nfig)
            it = 1:self.n_it;
            figure(nfig);
            hold on;
            plot(it, self.changes, "o-"); 
            title("CHANGES");
            hold off;
        end
        %-------------------------------------------
    end
    % END PUBLIC METHODS
    %----------------------------------------------------------------------  
      
    %----------------------------------------------------------------------
    %% PROTECTED METHODS 
    %----------------------------------------------------------------------
    methods (Access = protected)
        %-------------------------------------------
        % FUNCTIONAL EVALUATION
        function self = eval(self)
            nameF = self.name + "/func.txt";
            name_no_w_F = self.name + "/no_weight_func.txt";
            nameC = self.name + "/changes.txt";
            nameV = self.name + "/valid.txt";

            self.functional = readMat(nameF);
            self.f0Init = self.functional(1,1);
            self.f0Init_out_box = self.functional(1,2);
            self.f0Init_parts = self.functional(1,3:end);
            self.functional = self.functional(2:end, :);
            self.func = self.functional(1:end,1);
            self.n_it = length(self.func);
            self.func_out_box = self.functional(1:end,2);
            self.func_alpha = self.functional(1:end,3);
            self.func_grad  = self.functional(1:end,4);
            self.func_vort  = self.functional(1:end,5);
            self.func_p     = self.functional(1:end,6);
            self.valid = discardVec(nameV);
            self.changes =  discardVec(nameC);
            self.no_w_functional = readMat(name_no_w_F);
            self.no_w_f_alpha = self.no_w_functional(:,1);
            self.no_w_f_grad = self.no_w_functional(:,2);
            self.no_w_f_vort = self.no_w_functional(:,3);
            self.no_w_f_p = self.no_w_functional(:,4);

            self.func_values(:,1) = self.func;
            self.func_values(:,2) = self.func_out_box;
            self.func_values(:,3) = self.func_alpha;
            self.func_values(:,4) = self.func_grad;
            self.func_values(:,5) = self.func_vort;
            self.func_values(:,6) = self.func_p;
        end
        %-------------------------------------------
    end
    % END PROTECTED METHODS
    %----------------------------------------------------------------------  
end































