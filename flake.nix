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
          python3
          python3Packages.pip
          python3Packages.virtualenv
          python3Packages.build
          python3Packages.twine
          stdenv.cc.cc.lib
      ];
      shellHook = ''
          export LD_LIBRARY_PATH="${pkgs.stdenv.cc.cc.lib}/lib"
          test -d .nix-venv || ${pkgs.python3.interpreter} -m venv .nix-venv
          source .nix-venv/bin/activate
          pip install matplotlib pandas ete4
          pip install --no-binary ":all:" numpy
        '';
    };
    explorationShell = pkgs.mkShell{
        buildInputs = with pkgs; [
            python3
            python3Packages.pip
            python3Packages.virtualenv
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
