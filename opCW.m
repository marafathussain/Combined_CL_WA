classdef opCW < opSpot

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = private)
       funHandle = []; % Multiplication function
    end % Properties

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Constructor
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function op = opCW(m,n,nbscales,nbangles,ttype,pat,tp)

          if nargin < 3, nbscales = max(1,ceil(log2(min(m,n)) - 3)); end;
          if nargin < 4, nbangles = 16;                              end;
          if nargin < 5, ttype = 'ME';                               end; 
          if nargin < 6, pat = 'p';                                  end;
          if nargin < 7, tp = 'directional';                         end;

          finest  = 1;
          is_real = 1;

          % Compute length of curvelet coefficient vector
          C1 = mefcv2(randn(m,n),m,n,nbscales,nbangles);

          hdr{1}{1} = size(C1{1}{1});
          cnc = prod(hdr{1}{1});
          for i = 2:nbscales
             nw = length(C1{i});
             hdr{i}{1} = size(C1{i}{1});
             hdr{i}{2} = size(C1{i}{nw/2+1});
             cnc = cnc + nw/2*prod(hdr{i}{1}) + nw/2*prod(hdr{i}{2});
          end
          
          CL = cnc;
          
          % Compute length of wavwatom coefficient vector
          C2 = mefwa2sym(randn(m,n),pat,tp);
          s = size(C2);
          cnw = length(C2(:));
          
          cn = cnc + cnw;
          
          % parameter
          parms = {m,n,cnc,cnw,CL,hdr,finest,nbscales,nbangles,is_real,ttype,pat,tp,s};
          fun   = @(x,mode) opCW_intrnl(parms{:},x,mode);

          % Construct operator
          op = op@opSpot('CW', cn, m*n);
          op.cflag     = ~is_real;
          op.funHandle = fun;
       end % Constructor

    end % Methods
       
 
    methods ( Access = protected )
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Multiply
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function y = multiply(op,x,mode)
          y = op.funHandle(x,mode);
       end % Multiply          

    end % Methods
   
end % Classdef


%=======================================================================


function y = opCW_intrnl(m,n,cnc,cnw,CL,hdr,ac,nbs,nba,is_real,ttype,pat,tp,s,x,mode)

    if mode == 1
       % Analysis mode
       %size(x)
       y1 = mefcv2(reshape(x,m,n),m,n,nbs,nba);
       y1 = spot.utils.fdct_c2v(y1,cnc);
       
       y2 = mefwa2sym(reshape(x,m,n),pat,tp);
       y2 = reshape(y2,cnw,1);
       
       y = [y1;y2];
       %size(y)
    else
       % Synthesis mode  
       %size(x)
       x1 = mefdct_v2c(x(1:CL,1),hdr,nba);
       y1 = meicv2(x1,m,n,nbs,nba);
       
       x2 = reshape(x(CL+1:end,1),s);
       y2 = meiwa2sym(x2,pat,tp);
       
       y = y1 + y2;
       
       if is_real
          y = real(y);
       end
       y = y(:);
    end
end
