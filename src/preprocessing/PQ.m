function img_PQ = PQ( img )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% SMPTE ST.2084 definition
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define parameters
pq_n  = 2610/16384; 
pq_m  = 2523/32; 
pq_c1 = 3424/4096; 
pq_c2 = 2413/128; 
pq_c3 = 2392/128;

img = img ./ 4000;

PQcurve = @(x)( ( ( pq_c1 + pq_c2.*(x.^pq_n) )./( 1 + pq_c3.*(x.^pq_n) ) ).^pq_m );

img_PQ = single( real( PQcurve(img) ) );