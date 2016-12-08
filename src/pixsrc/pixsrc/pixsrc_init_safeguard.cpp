// check for possible setting configuration user did not really want
{
    double dummy;
    bool chk, thiskillpixsrc = 0;
    string msg = "Possible unintentional setting: ";
    
    // check regularization order
    chk = (-1==data_[g].regorder && data_[g].use_shapelets) ||
	0==data_[g].regorder || 1==data_[g].regorder || 2==data_[g].regorder;
    if (!chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "regularization order", cdata_->print2screenmutex);	
    }
    // check NPL stepsize
    chk = cdata_->npl_stepsize>=0;
    if (!chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "non-parametric lens step-size", cdata_->print2screenmutex);
    }
    // check NPL tolerance
    chk = cdata_->npl_ftolsize>0;
    if (!chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "non-parametric lens ftol (convergence criterion)", cdata_->print2screenmutex);
    }
    // check NPL TPS regularization stepsize
    chk = cdata_->npl_reg_parms[1]>=0;
    if (!chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "non-parametric lens TPS regularization stepsize", cdata_->print2screenmutex);
    }
    // check parametric source parameter estimation
    chk = 0==data_[g].guess_src_position || 1==data_[g].guess_src_position;
    if (!chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "analytic source parameter estimation", cdata_->print2screenmutex);
    }
    // source regularization flag
    chk = 0==data_[g].reg || 1==data_[g].reg;
    if (!chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "source regularization flag", cdata_->print2screenmutex);
    }
    // analytic source regularization
    chk = data_[g].usersetsrc && data_[g].reg;
    if (chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "you are using source regularization while fitting a parametric source to the data. "
			      "The gridded source will be regularized by the analytic source. Did you intend to do this?", cdata_->print2screenmutex);
    }
    // analytic source regularization with shapelets
    chk = data_[g].usersetsrc && data_[g].reg && data_[g].use_shapelets;
    if (chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "you are using analytic source regularization with shapelets turned on. "
			      "This is not supported and will lead to undefined behavior or a segmentation fault.", cdata_->print2screenmutex);
    }
    // srcexact with regularization
    chk = data_[g].usersetsrc && data_[g].reg && data_[g].srcexact;
    if (chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "you are using analytic source regularization with exact lens modeling of source. "
			      "This should not be done and will lead to undefined behavior or a segmentation fault.", cdata_->print2screenmutex);
    }
    // analytic source (without exact model) needs grid
    chk = data_[g].usersetsrc && !data_[g].srcexact && data_[g].use_shapelets;
    if (chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "analytic source (without exact modeling) needs grid. "
			      "Turn shapelets off, or this will lead to undefined behavior or a segmentation fault.", cdata_->print2screenmutex);
    }
    // source regularization and statistic mismatch
    chk = data_[g].reg && (data_[g].noevinochi || data_[g].onlychi);
    if (chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "you are using source regularization but not computing Bayesian evidence. "
			      "Did you intend to do this?", cdata_->print2screenmutex);
    }
    // source regularization and statistic mismatch
    chk = !data_[g].reg && (data_[g].noevinochi || !data_[g].onlychi);
    if (chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "you are not using source regularization and not computing chi^2 statistic. "
			      "Did you intend to do this?", cdata_->print2screenmutex);
    }
    // source regularization and statistic mismatch
    chk = data_[g].reg && data_[g].onlychi;
    if (chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "you are using source regularization, but you are computing chi^2 statistic. "
			      "Did you intend to do this?", cdata_->print2screenmutex);
    }
    // source regularization and statistic mismatch
    chk = !data_[g].reg && (!data_[g].noevinochi && !data_[g].onlychi);
    if (chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "you are not using source regularization, but you are computing Bayesian evidence. "
			      "Did you intend to do this?", cdata_->print2screenmutex);
    }
    // Ray direction (for creating lensing operator)
    chk = 1==data_[g].raydirection;
    if (!chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "ray direction should not be changed", cdata_->print2screenmutex);
    }
    // Numerical precision for ASCII files
    chk = data_[g].precision>0;
    if (!chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "numerical precision (ASCII files) should be positive", cdata_->print2screenmutex);
    }
    // magnification with shapelets
    chk = !data_[g].magparams || (data_[g].magparams && !data_[g].use_shapelets);
    if (!chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "as coded, computing magnification with shapelets can be dangerous because of the infinite extent of shapelets in source plane", cdata_->print2screenmutex);
    }
    // full magnification
    chk = 0==data_[g].fullmag || 1==data_[g].fullmag;
    if (!chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "full magnification option", cdata_->print2screenmutex);
    }
    // source pixel scale
    chk = !data_[g].srcpixelscale || data_[g].srcpixelscale[0]>0;
    if (!chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "source pixel scale must be positive", cdata_->print2screenmutex);
    }
    // source reg ftol
    chk = data_[g].regftol>0;
    if (!chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "tolerance for source regularization optimization must be positive", cdata_->print2screenmutex);
    }
    // shapelets
    chk = 0==data_[g].use_shapelets || 1==data_[g].use_shapelets;
    if (!chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "use shapelets option", cdata_->print2screenmutex);
    }
    // number of shapelets 
    chk = !data_[g].use_shapelets || (data_[g].num_shapelets1>0 && data_[g].num_shapelets2>0);
    if (!chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "number of shapelets must be positive", cdata_->print2screenmutex);
    }
    // number of shapelets 
    chk = !data_[g].use_shapelets || (data_[g].num_shapelets1==data_[g].num_shapelets2);
    if (!chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "problems could arise if number of shapelets in x and y direction are unequal", cdata_->print2screenmutex);
    }
    // shapelets scale
    chk = !data_[g].fixedshapeletparms[0] || data_[g].fixedshapeletparms[3]>0;
    if (!chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "shapelets scale must be positive", cdata_->print2screenmutex);
    }
    // shapelets pixel splitting
    chk = data_[g].shapeletpixsplit[0]>0 && data_[g].shapeletpixsplit[1]>=data_[g].shapeletpixsplit[0];
    if (!chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "shapelets pixel splitting must be positive and upper limit must be greater than or equal to lower limit", cdata_->print2screenmutex);
    }
    // source restart
    chk = data_[g].srcrestart>=1;
    if (!chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "source restart must be >= 1", cdata_->print2screenmutex);
    }
    // remove pixels
    chk = 0==data_[g].rmpix || 1==data_[g].rmpix;
    if (!chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "remove pixels option", cdata_->print2screenmutex);
    }
    // interperr
    chk = data_[g].interperr>=0;
    if (!chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "interp error parameter must be non-negative", cdata_->print2screenmutex);
    }
    // irscheme
    chk = 0==data_[g].irscheme || 1==data_[g].irscheme;
    if (!chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "interp error scheme parameter must be 0 or 1", cdata_->print2screenmutex);
    }
    // source ftol
    chk = data_[g].srcftol>0;
    if (!chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "source ftol must be positive", cdata_->print2screenmutex);
    }
    // source exact
    chk = 0==data_[g].srcexact || 1==data_[g].srcexact;
    if (!chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "exact source modeling parameter", cdata_->print2screenmutex);
    }
    // verbosity
    chk = data_[g].verbose>0;
    if (!chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "verbosity parameter must be positive", cdata_->print2screenmutex);
    }
    // source grid
    chk = 0==data_[g].gridtype || 1==data_[g].gridtype || 2==data_[g].gridtype;
    if (!chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "grid type must be one of 0,1,2", cdata_->print2screenmutex);
    }
    // source grid
    chk = 0==data_[g].gridtype;
    if (chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "currently, Cartesian grid is disabled", cdata_->print2screenmutex);
    }
    // source grid
    chk = 1==data_[g].gridtype && (data_[g].levelshift<0.5 || !OPERA equalszero(std::modf(data_[g].levelshift,&dummy)));
    if (chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "skip level of adaptive grid must be positive integer", cdata_->print2screenmutex);
    }
    // source grid
    chk = 2==data_[g].gridtype && data_[g].levelshift<=1.01;
    if (chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "irregular Cartesian grid shift level must be >1.01", cdata_->print2screenmutex);
    }
    // PSF
    chk = !data_[g].nopsf && !data_[g].psffromfile && !data_[g].majoraxis>=data_[g].minoraxis;
    if (chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "PSF major axis must be greater than or equal to minor axis", cdata_->print2screenmutex);
    }
    // PSF sampling
    chk = !data_[g].nopsf && !data_[g].psffromfile && !data_[g].psf_oversample>=1;
    if (chk) {
	thiskillpixsrc=1;
	PRINTER printwarning ("pixsrc", msg + "PSF oversampling must be >=1", cdata_->print2screenmutex);
    }
    
    
    killpixsrc = killpixsrc || (thiskillpixsrc && data_[g].fatalwarn);
}
