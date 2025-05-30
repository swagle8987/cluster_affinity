{
  description = "A very basic flake";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = inputs@{ self, nixpkgs, flake-utils }:
  flake-utils.lib.eachDefaultSystem (system:
  let 
    pkgs = import nixpkgs {inherit system;};
  in
    {
    devShell = pkgs.mkShell {
      buildInputs = with pkgs; [
          python310
          python310Packages.pip
          python310Packages.virtualenv
          python310Packages.build
          python310Packages.twine
          stdenv.cc.cc.lib
      ];
      shellHook = ''
          export LD_LIBRARY_PATH="${pkgs.stdenv.cc.cc.lib}/lib"
          test -d .nix-venv || ${pkgs.python3.interpreter} -m venv .nix-venv
          source .nix-venv/bin/activate
          pip install dendropy matplotlib pandas ete4
          pip install --no-binary ":all:" numpy
        '';
    };
    explorationShell = pkgs.mkShell{
        buildInputs = with pkgs; [
            python310
            python310Packages.pip
            python310Packages.virtualenv
            jupyter-all
            cmake
      ];
      shellHook = ''
          test -d .nix-venv || ${pkgs.python3.interpreter} -m venv .nix-venv
          source .nix-venv/bin/activate
        '';
    };
  });
}
