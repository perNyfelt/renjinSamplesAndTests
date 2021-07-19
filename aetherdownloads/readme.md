# Aetherdownloads

This example shows a bit of the details how the AetherPackageLoader works.

The AetherPackageLoader is a Renjin feature that omits the need to have all packages defined in the pom.xml
or equivalent. Very convenient sometimes but on the other hand there is much less control over what is actually 
running. See [package loading](https://renjin.readthedocs.io/en/latest/library/execution-context.html?highlight=AetherPackageLoader#package-loading)
in the renjin documentation for more info.