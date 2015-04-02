from _gfs_dynamics import _gfs_dynamics

test = _gfs_dynamics(384,190);

test.initModel();

test.oneStepForward();

test.configureModel(192,94);

test.initModel();

test.oneStepForward();

test.shutDownModel();
