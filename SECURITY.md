# Security policy

## Reporting a vulnerability

Please report it **privately**, not in a public issue.

Go to the **Security** tab of this repository, choose **Report a vulnerability**, and fill in the form. Only you and the maintainer can see the report.

Please include:

- what happens and why you think it is a security problem,
- the steps to reproduce it (the entry and its inputs, or the command line and the file that triggers it),
- the plugin version (the zip's name), your SolidCraft version and your operating system.

You will get an acknowledgement as soon as I see the report. For a verified problem I aim to have a fix released within about 90 days, and I will tell you if it will take longer. If you want to be credited in the release notes, say so.

## Supported versions

Only the latest release gets security fixes. Please check that the problem is still there in the latest release before you report it.

## What is in scope

- the plugin (`src/`, `wit/`, `plugin.json`): anything an input can make it do beyond what its manifest grants (it asks for the `commands@0.1` capability only), or that crashes or hangs SolidCraft;
- the command-line tools (`stress3d`, `stress3d_printed`, `beam_opt`): memory-safety or resource-exhaustion problems from a crafted STL file or arguments;
- the build and CI scripts and the files attached to a release (the zip and `SHA256SUMS.txt`).

## What is not a vulnerability by itself

- A problem in SolidCraft itself belongs to [tangency-labs/solidcraft](https://github.com/tangency-labs/solidcraft).
- A wrong number (a material value, a stress, a safety factor) is a bug, not a vulnerability: open a normal issue. The results are a design aid, not a certified analysis; see the README.

## Disclosure

Please give me the time above before you publish details. When a fix is released I will publish a GitHub security advisory that names the affected versions.
