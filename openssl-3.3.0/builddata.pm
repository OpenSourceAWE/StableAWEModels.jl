package OpenSSL::safe::installdata;

use strict;
use warnings;
use Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw($PREFIX
                  $BINDIR $BINDIR_REL
                  $LIBDIR $LIBDIR_REL
                  $INCLUDEDIR $INCLUDEDIR_REL
                  $APPLINKDIR $APPLINKDIR_REL
                  $ENGINESDIR $ENGINESDIR_REL
                  $MODULESDIR $MODULESDIR_REL
                  $PKGCONFIGDIR $PKGCONFIGDIR_REL
                  $CMAKECONFIGDIR $CMAKECONFIGDIR_REL
                  $VERSION @LDLIBS);

our $PREFIX             = '/home/jellepoland/ownCloud/phd/code/OpenSourceAWE/SymbolicAWEModels.jl/openssl-3.3.0';
our $BINDIR             = '/home/jellepoland/ownCloud/phd/code/OpenSourceAWE/SymbolicAWEModels.jl/openssl-3.3.0/apps';
our $BINDIR_REL         = 'apps';
our $LIBDIR             = '/home/jellepoland/ownCloud/phd/code/OpenSourceAWE/SymbolicAWEModels.jl/openssl-3.3.0';
our $LIBDIR_REL         = '.';
our $INCLUDEDIR         = '/home/jellepoland/ownCloud/phd/code/OpenSourceAWE/SymbolicAWEModels.jl/openssl-3.3.0/include';
our $INCLUDEDIR_REL     = 'include';
our $APPLINKDIR         = '/home/jellepoland/ownCloud/phd/code/OpenSourceAWE/SymbolicAWEModels.jl/openssl-3.3.0/ms';
our $APPLINKDIR_REL     = 'ms';
our $ENGINESDIR         = '/home/jellepoland/ownCloud/phd/code/OpenSourceAWE/SymbolicAWEModels.jl/openssl-3.3.0/engines';
our $ENGINESDIR_REL     = 'engines';
our $MODULESDIR         = '/home/jellepoland/ownCloud/phd/code/OpenSourceAWE/SymbolicAWEModels.jl/openssl-3.3.0/providers';
our $MODULESDIR_REL     = 'providers';
our $PKGCONFIGDIR       = '';
our $PKGCONFIGDIR_REL   = '';
our $CMAKECONFIGDIR     = '';
our $CMAKECONFIGDIR_REL = '';
our $VERSION            = '3.3.0';
our @LDLIBS             =
    # Unix and Windows use space separation, VMS uses comma separation
    split(/ +| *, */, '-ldl -pthread ');

1;
