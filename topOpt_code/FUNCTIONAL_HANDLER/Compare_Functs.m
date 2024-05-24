
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
            for ifunc=1:self.n_func
                functional = self.functionals(ifunc);
                func = functional.func;
                it = 1:functional.n_it;
                plot(it, func, "-","LineWidth", 2, 'DisplayName', functional.custom_name);
            end
            legend;
            colororder("reef");
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
