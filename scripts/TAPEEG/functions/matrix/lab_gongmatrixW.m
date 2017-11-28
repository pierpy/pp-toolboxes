% Output is weighted connectivitiy matrix of 78 source regions according to
% Gong et al.

function matrix = lab_gongmatrixW

matrix = [0,0.842200000000000,0.934000000000000,0.863500000000000,0,0.770400000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.805500000000000,0,0,0,0.923000000000000,0,0,0.836700000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0.842200000000000,0,0.814000000000000,0.767700000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.814900000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0.934000000000000,0.814000000000000,0,0.883700000000000,0.834600000000000,0.823800000000000,0.601500000000000,0,0,0,0.720600000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.692300000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0.863500000000000,0.767700000000000,0.883700000000000,0,0,0,0,0,0,0,0.817600000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.870300000000000,0,0,0,0,0,0,0.903800000000000,0,0,0,0,0,0,0.791600000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0.834600000000000,0,0,0.881400000000000,0.630800000000000,0.780600000000000,0,0.801000000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0.770400000000000,0,0.823800000000000,0,0.881400000000000,0,0,0,0,0.840200000000000,0,0,0,0,0,0,0,0,0,0,0,0,0.230900000000000,0.242100000000000,0,0,0,0.414300000000000,0,0,0,0,0.868200000000000,0,0,0,0,0,0.813500000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0.601500000000000,0,0.630800000000000,0,0,0.816600000000000,0,0.718600000000000,0.789900000000000,0.739600000000000,0.487600000000000,0.630400000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.762700000000000,0,0,0,0.747200000000000,0.701000000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0.780600000000000,0,0.816600000000000,0,0.765100000000000,0.876800000000000,0.763400000000000,0.602000000000000,0,0.564800000000000,0,0.618900000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0.765100000000000,0,0.848500000000000,0,0,0,0.443600000000000,0.626600000000000,0.567800000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.527000000000000,0,0,0,0,0,0,0,0.822100000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0.801000000000000,0.840200000000000,0.718600000000000,0.876800000000000,0.848500000000000,0,0,0,0,0.563500000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.833400000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0.720600000000000,0.817600000000000,0,0,0.789900000000000,0.763400000000000,0,0,0,0.559800000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.871700000000000,0.549600000000000,0,0,0,0,0,0.802900000000000,0,0,0.747200000000000,0,0,0,0.916600000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.843400000000000,0.544200000000000,0,0;0,0,0,0,0,0,0.739600000000000,0.602000000000000,0,0,0.559800000000000,0,0.738300000000000,0.800700000000000,0,0,0,0,0,0,0.671300000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.850300000000000,0,0,0,0,0,0,0,0,0,0,0,0,0.549600000000000,0.903800000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0.487600000000000,0,0,0,0,0.738300000000000,0,0.777500000000000,0,0.616400000000000,0,0,0,0,0.898400000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.745700000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0.719000000000000,0.910200000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0.630400000000000,0.564800000000000,0.443600000000000,0.563500000000000,0,0.800700000000000,0.777500000000000,0,0.701200000000000,0.794700000000000,0.716300000000000,0.747700000000000,0.697500000000000,0.636100000000000,0,0,0,0,0,0,0,0,0,0,0.507400000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0.626600000000000,0,0,0,0,0.701200000000000,0,0.855700000000000,0,0,0.825800000000000,0,0,0,0,0,0,0,0,0,0.941600000000000,0.880700000000000,0,0,0,0,0,0,0,0,0.803800000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0.618900000000000,0.567800000000000,0,0,0,0.616400000000000,0.794700000000000,0.855700000000000,0,0.604000000000000,0.777600000000000,0.831500000000000,0.694200000000000,0.628200000000000,0,0,0,0,0,0,0,0,0.795000000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0.716300000000000,0,0.604000000000000,0,0.770400000000000,0,0.717100000000000,0.891700000000000,0.708400000000000,0.617600000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.800300000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0.747700000000000,0,0.777600000000000,0.770400000000000,0,0.870100000000000,0.881700000000000,0,0,0.704300000000000,0,0,0,0,0,0,0,0.661300000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0.697500000000000,0.825800000000000,0.831500000000000,0,0.870100000000000,0,0.841100000000000,0,0,0,0,0,0,0,0,0,0.799200000000000,0.774100000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0.636100000000000,0,0.694200000000000,0.717100000000000,0.881700000000000,0.841100000000000,0,0,0,0.817600000000000,0,0,0,0,0,0,0,0.702400000000000,0.619400000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0.671300000000000,0.898400000000000,0,0,0.628200000000000,0.891700000000000,0,0,0,0,0.711300000000000,0.623300000000000,0,0.664200000000000,0.712700000000000,0.572200000000000,0,0,0,0,0,0,0,0.368600000000000,0.338600000000000,0.725300000000000,0.771400000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0.852200000000000,0,0,0,0,0,0,0,0.884600000000000,0.642800000000000,0,0,0.630800000000000,0.683300000000000,0,0,0,0,0,0,0,0,0,0,0,0.746800000000000,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.708400000000000,0,0,0,0.711300000000000,0,0.854300000000000,0,0.821600000000000,0.914400000000000,0.713100000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.668300000000000,0.724200000000000,0,0,0.726800000000000,0.804000000000000,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0.230900000000000,0,0,0,0,0,0,0,0,0,0,0.617600000000000,0.704300000000000,0,0.817600000000000,0.623300000000000,0.854300000000000,0,0.819800000000000,0.827700000000000,0.806900000000000,0.778800000000000,0.734200000000000,0,0.557100000000000,0.643500000000000,0.599600000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.607700000000000,0,0,0.673300000000000,0.696600000000000,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0.242100000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.819800000000000,0,0.785300000000000,0,0.833700000000000,0.823800000000000,0,0,0.659600000000000,0.665400000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.664200000000000,0.821600000000000,0.827700000000000,0.785300000000000,0,0.861000000000000,0.878100000000000,0,0,0,0,0,0,0,0.499000000000000,0,0,0.809300000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.640700000000000,0.711700000000000,0,0,0.833200000000000,0.799800000000000,0.791800000000000,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.712700000000000,0.914400000000000,0.806900000000000,0,0.861000000000000,0,0.740700000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.697700000000000,0.797700000000000,0,0,0.809100000000000,0.884600000000000,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.572200000000000,0.713100000000000,0.778800000000000,0.833700000000000,0.878100000000000,0.740700000000000,0,0.879200000000000,0,0,0,0.689400000000000,0,0,0.585300000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0.414300000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.734200000000000,0.823800000000000,0,0,0.879200000000000,0,0,0,0.778800000000000,0.804800000000000,0.522600000000000,0.512800000000000,0.637000000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.941600000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.909100000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.880700000000000,0.795000000000000,0,0,0.799200000000000,0,0,0,0.557100000000000,0,0,0,0,0,0.909100000000000,0,0.853000000000000,0,0.720300000000000,0,0,0,0,0,0.798200000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0.527000000000000,0,0,0,0,0.507400000000000,0,0,0,0.661300000000000,0.774100000000000,0.702400000000000,0,0,0.643500000000000,0.659600000000000,0,0,0,0.778800000000000,0,0.853000000000000,0,0.894200000000000,0.632600000000000,0.590300000000000,0,0,0,0,0.678500000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.619400000000000,0,0,0.599600000000000,0.665400000000000,0,0,0.689400000000000,0.804800000000000,0,0,0.894200000000000,0,0,0.632600000000000,0.677200000000000,0,0,0,0.659100000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0.868200000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.522600000000000,0,0.720300000000000,0.632600000000000,0,0,0.898600000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.512800000000000,0,0,0.590300000000000,0.632600000000000,0.898600000000000,0,0.807600000000000,0,0,0,0.769600000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.368600000000000,0,0,0,0.499000000000000,0,0.585300000000000,0.637000000000000,0,0,0,0.677200000000000,0,0.807600000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0.805500000000000,0.814900000000000,0,0.870300000000000,0,0,0,0,0,0,0.871700000000000,0,0,0,0,0,0,0,0,0,0.338600000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.610300000000000,0.433800000000000,0,0,0,0,0,0,0,0,0,0,0,0.847000000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.910200000000000,0,0,0;0,0,0,0,0,0,0,0,0,0,0.549600000000000,0.850300000000000,0.745700000000000,0,0,0,0,0,0,0,0.725300000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.610300000000000,0,0.741800000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.771400000000000,0,0,0,0.809300000000000,0,0,0,0,0,0,0,0,0,0,0.433800000000000,0.741800000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.761900000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.923000000000000,0;0,0,0.692300000000000,0,0,0.813500000000000,0,0,0.822100000000000,0.833400000000000,0,0,0,0,0.803800000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.798200000000000,0.678500000000000,0.659100000000000,0,0.769600000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0.923000000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.842900000000000,0.922500000000000,0.863000000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.842900000000000,0,0.811300000000000,0.767700000000000,0,0,0,0,0,0,0.685200000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.704100000000000,0,0.811000000000000,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.922500000000000,0.811300000000000,0,0.880700000000000,0.840000000000000,0.835400000000000,0.600700000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.745100000000000,0,0,0,0,0;0.836700000000000,0,0,0.903800000000000,0,0,0,0,0,0,0.802900000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.863000000000000,0.767700000000000,0.880700000000000,0,0,0,0,0,0,0,0.818500000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.871100000000000,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.840000000000000,0,0,0.882300000000000,0,0.780600000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.241300000000000,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.835400000000000,0,0.882300000000000,0,0,0,0,0.834500000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.350200000000000,0,0,0,0,0,0,0,0,0,0,0,0.811800000000000;0,0,0,0,0,0,0.762700000000000,0,0,0,0.747200000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.600700000000000,0,0,0,0,0.834500000000000,0.628700000000000,0.727200000000000,0.789900000000000,0.736100000000000,0,0.630400000000000,0,0.588900000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.763300000000000,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.780600000000000,0,0.834500000000000,0,0.765700000000000,0.870300000000000,0.786000000000000,0.608200000000000,0,0.565900000000000,0,0.609600000000000,0,0.412600000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.628700000000000,0.765700000000000,0,0.843300000000000,0,0,0,0.446000000000000,0.623300000000000,0.559900000000000,0,0,0.465400000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.821100000000000;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.834500000000000,0.727200000000000,0.870300000000000,0.843300000000000,0,0,0,0,0.565900000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.833400000000000;0,0,0,0.791600000000000,0,0,0.747200000000000,0,0,0,0.916600000000000,0.549600000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.847000000000000,0,0,0,0,0.685200000000000,0,0.818500000000000,0,0,0.789900000000000,0.786000000000000,0,0,0,0.560000000000000,0,0,0,0,0,0,0,0,0.277800000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.871600000000000,0.549500000000000,0,0;0,0,0,0,0,0,0.701000000000000,0,0,0,0,0.903800000000000,0.719000000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.736100000000000,0.608200000000000,0,0,0.560000000000000,0,0.738100000000000,0.789800000000000,0,0,0,0,0,0,0.672300000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.849900000000000,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0.910200000000000,0,0,0,0,0,0,0,0.852200000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.738100000000000,0,0.753400000000000,0,0.586600000000000,0,0,0,0,0.900000000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.746000000000000,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.630400000000000,0.565900000000000,0.446000000000000,0.565900000000000,0,0.789800000000000,0.753400000000000,0,0.699600000000000,0.790600000000000,0.710500000000000,0.758800000000000,0.712500000000000,0.637500000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.623300000000000,0,0,0,0,0.699600000000000,0,0.855700000000000,0,0,0.840600000000000,0,0,0,0,0,0,0,0,0,0,0.900200000000000,0,0,0,0,0,0,0,0,0.801600000000000;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.588900000000000,0.609600000000000,0.559900000000000,0,0,0,0.586600000000000,0.790600000000000,0.855700000000000,0,0.590100000000000,0.761900000000000,0.833000000000000,0.692000000000000,0.589500000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.710500000000000,0,0.590100000000000,0,0.787500000000000,0.650700000000000,0.709600000000000,0.874800000000000,0.708200000000000,0.615900000000000,0,0,0,0,0,0,0,0.455400000000000,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.412600000000000,0,0,0,0,0,0.758800000000000,0,0.761900000000000,0.787500000000000,0,0.862700000000000,0.875600000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.465400000000000,0,0,0,0,0.712500000000000,0.840600000000000,0.833000000000000,0.650700000000000,0.862700000000000,0,0.850700000000000,0,0,0,0,0,0,0,0,0,0.798000000000000,0.774100000000000,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.637500000000000,0,0.692000000000000,0.709600000000000,0.875600000000000,0.850700000000000,0,0,0.777600000000000,0.815200000000000,0,0,0,0,0,0,0.672100000000000,0.709700000000000,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.800300000000000,0,0,0,0.884600000000000,0.668300000000000,0,0,0.640700000000000,0.697700000000000,0,0,0,0,0,0,0,0,0,0,0,0.761900000000000,0,0,0,0,0,0,0,0,0,0,0,0.277800000000000,0.672300000000000,0.900000000000000,0,0,0.589500000000000,0.874800000000000,0,0,0,0,0.699100000000000,0.608600000000000,0,0.662200000000000,0.712700000000000,0.572900000000000,0,0,0,0,0,0,0,0.364500000000000,0.339500000000000,0.727100000000000,0.774700000000000,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.642800000000000,0.724200000000000,0.607700000000000,0,0.711700000000000,0.797700000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.708200000000000,0,0,0.777600000000000,0.699100000000000,0,0.863000000000000,0.703900000000000,0.815800000000000,0.908600000000000,0.707400000000000,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.615900000000000,0,0,0.815200000000000,0.608600000000000,0.863000000000000,0,0.820000000000000,0.822100000000000,0.812900000000000,0.757500000000000,0.732500000000000,0,0.559200000000000,0.651700000000000,0.603700000000000,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.703900000000000,0.820000000000000,0,0.776600000000000,0,0.802800000000000,0.821000000000000,0,0,0,0.673100000000000,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.630800000000000,0.726800000000000,0.673300000000000,0,0.833200000000000,0.809100000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.662200000000000,0.815800000000000,0.822100000000000,0.776600000000000,0,0.863300000000000,0.884200000000000,0.777300000000000,0,0,0,0,0,0,0.498300000000000,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.683300000000000,0.804000000000000,0.696600000000000,0,0.799800000000000,0.884600000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.712700000000000,0.908600000000000,0.812900000000000,0,0.863300000000000,0,0.749100000000000,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.791800000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.241300000000000,0.350200000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.572900000000000,0.707400000000000,0.757500000000000,0.802800000000000,0.884200000000000,0.749100000000000,0,0.863300000000000,0,0,0,0,0,0,0.582400000000000,0,0,0.769600000000000,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.732500000000000,0.821000000000000,0.777300000000000,0,0.863300000000000,0,0,0,0.778800000000000,0.799900000000000,0,0.512800000000000,0.640900000000000,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.909100000000000,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.900200000000000,0,0,0,0.798000000000000,0.672100000000000,0,0,0.559200000000000,0,0,0,0,0,0.909100000000000,0,0.853600000000000,0.783600000000000,0.723200000000000,0.644100000000000,0,0,0,0,0.801600000000000;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.455400000000000,0,0.774100000000000,0.709700000000000,0,0,0.651700000000000,0,0,0,0,0.778800000000000,0,0.853600000000000,0,0.895600000000000,0.637600000000000,0.590300000000000,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.603700000000000,0.673100000000000,0,0,0,0.799900000000000,0,0.783600000000000,0.895600000000000,0,0,0.631000000000000,0.686300000000000,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.723200000000000,0.637600000000000,0,0,0.887800000000000,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.704100000000000,0.745100000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.512800000000000,0,0.644100000000000,0.590300000000000,0.631000000000000,0.887800000000000,0,0.821200000000000,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.364500000000000,0,0,0,0.498300000000000,0,0.582400000000000,0.640900000000000,0,0,0,0.686300000000000,0,0.821200000000000,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0.843400000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.910200000000000,0,0,0,0,0.811000000000000,0,0.871100000000000,0,0,0.763300000000000,0,0,0,0.871600000000000,0,0,0,0,0,0,0,0,0,0.339500000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.609900000000000,0,0;0,0,0,0,0,0,0,0,0,0,0.544200000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.549500000000000,0.849900000000000,0.746000000000000,0,0,0,0,0,0,0,0.727100000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.609900000000000,0,0.741800000000000,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.746800000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.923000000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.774700000000000,0,0,0,0,0,0.769600000000000,0,0,0,0,0,0,0,0,0,0.741800000000000,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.811800000000000,0,0,0.821100000000000,0.833400000000000,0,0,0,0,0.801600000000000,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.801600000000000,0,0,0,0,0,0,0,0,0];