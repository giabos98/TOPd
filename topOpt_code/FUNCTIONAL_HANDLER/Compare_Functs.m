
%--------------------------------------------------------------------------
% CLASS CREATED TO COMPARE FUNCTIONALS
%--------------------------------------------------------------------------

classdef Compare_Functs
    %----------------------------------------------------------------------
    %% PROPERTIES
    %----------------------------------------------------------------------
    properties
        % general properties
        n_func = 0;
        functionals;
        min_it = 1e10;
    end
    % END PROPERTIES
    %----------------------------------------------------------------------

    %----------------------------------------------------------------------
    %% PUBLIC METHODS 
    %----------------------------------------------------------------------
    methods (Access = public)
        %-------------------------------------------
        % CONSTRUCTOR
        function self = Compare_Functs() 
            self.n_func = 0;
            self.functionals = [];
        end
        %-------------------------------------------

        %-------------------------------------------
        % ADD FUNCTIONAL
        function self = add_functional(self, functional)
            self.functionals = [self.functionals, functional];
            self.n_func = self.n_func+1;
            if (functional.n_it < self.min_it)
                self.min_it = functional.n_it;
            end
        end
        %-------------------------------------------

        %-------------------------------------------
        % PRINT FUNCTIONALS
        function new_nfig = print(self, nfig, what_to_print)
            % print starting from nfig+1
            % what_to_print = [print_func, print_changes, print_no_w_func]
            new_nfig = nfig;
            for ifunc=1:self.n_func
                new_nfig = self.functionals(ifunc).print(new_nfig, what_to_print);
            end
        end
        %-------------------------------------------

        %-------------------------------------------
        % COMPARE FUNCTIONALS
        function new_nfig = compare_base_functionals(self, nfig, normalize, print_changes, axis_scale)
            % print in nfig+1
            new_nfig = nfig+1;
            figure(new_nfig);
            scales = ones(self.n_func,1);
            % normalize = 1: use functional values relative to their first
            % one
            if (normalize(1) ~= 1)
                for ifunc=1:self.n_func
                    scales(ifunc) = self.functionals(ifunc).f0Init;
                end
                if (normalize(1) == -1) 
                    if (normalize(2) < 1e-10)
                        rescale = 0;
                        %set rescale to max(f0Init)
                        for ifunc=1:self.n_func
                            if (self.functionals(ifunc).f0Init > rescale)
                                rescale = self.functionals(ifunc).f0Init;
                            end
                        end
                    else %rescale for a given factor > 0
                        rescale = normalize(2);
                    end
                    scales = scales / rescale;
                elseif (normalize(1) == -2)
                    for ifunc=1:self.n_func
                        scales(ifunc) =  scales(ifunc) / self.functionals(ifunc).mu;
                    end
                end
            end
            
            hold on;
            title("BASE FUNCTIONAL COMPARISON");
            it = 1:self.min_it;
            for ifunc=1:self.n_func
                functional = self.functionals(ifunc);
                func = functional.func;
                func = [1;func(1:self.min_it)];
                func = func * scales(ifunc);
                plot([0,it], func, "-","LineWidth", 2, 'DisplayName', functional.custom_name);
                colororder("reef");
            end

            if (print_changes == 1)
                for ifunc=1:self.n_func
                    functional = self.functionals(ifunc);
                    func = functional.func;
                    func = func * scales(ifunc);
                    valid = functional.valid;
                    for i = 1 : floor(self.min_it / self.min_it) : (self.min_it)
                        if valid(i) == 1
                            plot(it(i), func(i), 'bo', "LineWidth", 0.5, 'HandleVisibility','off');
                        else
                            plot(it(i), func(i), 'ro', "LineWidth", 0.5, 'HandleVisibility','off');
                        end
                    end
                end   
            end
            if (axis_scale == "log")
                xscale("log");
            end
            legend;
            
            hold off;
        end
        %-------------------------------------------
        % COMPARE FUNCTIONALS
        function new_nfig = compare_functionals(self, nfig, value_to_compare, abs_value, print_changes, axis_scale, custom_colors, max_it)
            % print in nfig+1
            % what_to_print = [print_func, print_changes, print_no_w_func]
            new_nfig = nfig+1;
            figure(new_nfig);
            
            if (length(custom_colors) == self.n_func)
                temp_colors = custom_colors;
            else
                for ifunc=1:self.n_func
                    temp_colors(ifunc) = self.get_random_color();
                end
            end
            
            hold on;
            % title("BASE FUNCTIONAL COMPARISON");
            for ifunc=1:self.n_func
                functional = self.functionals(ifunc);
                functional.print_specific_func(new_nfig, value_to_compare, abs_value, print_changes, temp_colors(ifunc), max_it);
            end

            if (axis_scale == "log")
                xscale("log");
            end
            % legend;
            
            hold off;
        end
        %---
        function new_nfig = compare_functionals_initial_values(self, nfig)
            new_nfig = nfig+1;
            x_label_loc = [];
            x_label = [];
            figure(new_nfig);
            hold on;
            for ifunc = 1:self.n_func
                temp_func = self.functionals(ifunc);
                rand_col = self.get_random_color();
                plot(ifunc, temp_func.f0Init, "o", "Color", rand_col);
                x_labels_loc = [x_label_loc, ifunc];
                x_labels = [x_label, temp_func.custom_name];
            end
            xticks(x_labels_loc)
            xticklabels(x_labels)

            hold off;
        end
        %-------------------------------------------

        function new_nfig = compare_alpha_max(self, nfig)
            new_nfig = nfig+1;
            for ifunc=1:self.n_func
                functional = self.functionals(ifunc);
                functional.print_alpha_max(new_nfig);
            end
        end
        %-------------------------------------------
    end
    % END PUBLIC METHODS
    %----------------------------------------------------------------------  
      
    %----------------------------------------------------------------------
    %% PROTECTED METHODS 
    %----------------------------------------------------------------------
    methods (Access = protected)
    end
    % END PROTECTED METHODS
    %----------------------------------------------------------------------  

    %----------------------------------------------------------------------
    %% STATIC METHODS 
    %----------------------------------------------------------------------
    methods (Static)
        function color = get_random_color()
            color = "#";
            r = dec2hex(randi([0,255]), 2);
            g = dec2hex(randi([0,255]), 2);
            b = dec2hex(randi([0,255]), 2);
            color = color + r + g + b;
        end
        %---

        function func_case = get_func_case(what_case)
            switch what_case
                case "tot"
                    func_case = 1;
                case "out"
                    func_case = 2;
                case "alpha"
                    func_case = 3;
                case "grad"
                    func_case = 4;
                case "vort"
                    func_case = 5;
                case "p"
                    func_case = 6;
                otherwise
                    func_case = 1;
            end
        end
        %---
    end
    % END STATIC METHODS
    %---------------------------------------------------------------------- 
end
