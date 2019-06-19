function imgOut = cat_0( imgOutR,   imgOutG, imgOutB )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Concatenates three channels discarding 0 values
%% Inputs: 1. imgOutR -> MxN red channel image (ídem for G and B)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

saturatedPixelsRed   = find( imgOutR == 0 );
saturatedPixelsGreen = find( imgOutG == 0 );
saturatedPixelsBlue  = find( imgOutB == 0 );

imgOutR( saturatedPixelsGreen ) = 0; imgOutR( saturatedPixelsBlue ) = 0;
imgOutG( saturatedPixelsRed )   = 0; imgOutG( saturatedPixelsBlue ) = 0;
imgOutB( saturatedPixelsGreen ) = 0; imgOutB( saturatedPixelsRed )  = 0;

imgOut = cat(3, imgOutR, imgOutG, imgOutB);