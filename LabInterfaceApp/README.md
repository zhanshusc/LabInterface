# Lab Interface application source

This directory is the clean source root for running and packaging the Python
version of the Lab Interface.

## Included

- `lab_interface.py`: application entry point
- `imports/TAExperiment.py`: transient-absorption processing
- `imports/SolitonTAExperimentSelfContained.py`: Soliton TA processing
- `requirements.txt`: pinned application dependencies
- `requirements-build.txt`: application dependencies plus PyInstaller
- `LabInterface-macos-arm64.spec`: native ARM64 PyInstaller configuration
- `build_macos_arm64.sh`: architecture checks, build, validation, and ZIP creation

The application opens experiment data through file dialogs, so example TA and
TCSPC datasets are intentionally not bundled with the application source.

## Excluded from this build source

- MATLAB source and MATLAB Engine integrations
- older Lab Interface implementations
- standalone TCSPC experiments
- test scripts and test datasets
- editor settings and generated Python cache files

Run the source version from this directory with:

```sh
python lab_interface.py
```

The PyInstaller configuration and GitHub Actions workflow use this directory
as their source root.

## Build through GitHub Actions

The repository workflow `.github/workflows/build-macos-arm64.yml` runs on an
Apple Silicon `macos-15` runner. It installs the pinned dependencies, builds a
native ARM64 application, verifies the executable architecture and ad-hoc code
signature, and uploads `LabInterface-macOS-arm64.zip` as a workflow artifact.

To build it:

1. Commit and push `LabInterfaceApp` and `.github/workflows` to GitHub.
2. Open the repository's **Actions** page.
3. Select **Build macOS Apple Silicon app**.
4. Choose **Run workflow**.
5. When the run finishes, download **LabInterface-macOS-arm64** from its
   **Artifacts** section.

This produces an ad-hoc-signed application suitable for testing. Public
distribution without macOS security warnings requires an Apple Developer ID
certificate and notarization, which are intentionally not configured here.
