
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
        function new_nfig = compare_base_functionals(self, nfig)
            % print in nfig+1
            % what_to_print = [print_func, print_changes, print_no_w_func]
            new_nfig = nfig+1;
            figure(new_nfig);
            hold on;
            title("BASE FUNCTIONAL COMPARISON");
            it = 1:self.min_it;
            for ifunc=1:self.n_func
                functional = self.functionals(ifunc);
                func = functional.func;
                func = func(1:self.min_it);
                plot(it, func, "-","LineWidth", 2, 'DisplayName', functional.custom_name);
                colororder("reef");
            end
            for ifunc=1:self.n_func
                functional = self.functionals(ifunc);
                func = functional.func;
                valid = functional.valid;
                for i = 1 : floor(self.min_it / 20) : (self.min_it)
                    if valid(i) == 1
                        plot(it(i), func(i), 'bo', "LineWidth", 1.5, 'HandleVisibility','off');
                    else
                        plot(it(i), func(i), 'ro', "LineWidth", 1.5, 'HandleVisibility','off');
                    end
                end
            end            
            legend;
            
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
    end
    % END PROTECTED METHODS
    %----------------------------------------------------------------------  
end
